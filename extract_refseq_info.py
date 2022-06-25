#!/usr/bin/env python

import json
import os
import sys
import traceback
from typing import Any, Dict, List, Union

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

FULL_ACCESSIONS = {}
with open("refseq_to_new_accessions") as handle:
    pairs = handle.read().splitlines()
    for pair in pairs:
        bgc, acc = pair.split()
        FULL_ACCESSIONS[bgc] = acc

EMPTIES = {}
with open("empty_refseq_dates") as handle:
    pairs = handle.read().splitlines()
    for pair in pairs:
        bgc, acc, date = pair.split()
        if os.sep in bgc:
            bgc = os.path.basename(bgc)
        if bgc.endswith(".json"):
            bgc.rsplit(".", 1)[0]
        EMPTIES[bgc] = (acc, pair)


def location_to_json(location: Union[CompoundLocation, FeatureLocation]) -> Dict[str, Any]:
    exons = [{"start": part.start, "end": part.end} for part in location.parts]
    return {"exons": exons, "strand": location.strand}


def get_refseq_base(comment: str) -> str:
    raw_sentences = [raw.strip() for raw in comment.replace("\n", " ").split(".")]
    sentences = raw_sentences[:1]
    for i, raw in enumerate(raw_sentences[1:]):
        if raw.isdigit():
            sentences[-1] = f"{sentences[-1]}.{raw}"
        else:
            sentences.append(raw)
    for sentence in sentences:
        if "reference sequence is identical to" in sentence or "reference sequence was derived from" in sentence:
            result = sentence.rsplit()[-1]
            return result
    raise ValueError(f"couldn't get refseq base accession from: {comment}")


def update_gene_naming(data: Dict[str, Any], name_mapping: Dict[str, str]) -> None:
    def replace(subdata, key):
        if key in subdata and subdata[key] in name_mapping:
            subdata[key] = name_mapping[subdata[key]]

    nrp = data["cluster"].get("nrp", {}) 
    if nrp:
        for nrps_gene in nrp.get("nrps_genes", []):
            replace(nrps_gene, "gene_id")
        for thioesterase in nrp.get("thioesterases", []):
            replace(thioesterase, "gene")
    polyketide = data["cluster"].get("polyketide", {})
    for synthase in polyketide.get("synthases", []):
        if "genes" not in synthase:
            continue
        synthase["genes"] = [name_mapping.get(gene, gene) for gene in synthase["genes"]]
        for module in synthase.get("modules", []):
            module["genes"] = [name_mapping.get(gene, gene) for gene in synthase["genes"]]
    cyclases = polyketide.get("cyclases", [])
    if cyclases:
         polyketide["cyclases"] = [name_mapping.get(gene, gene) for gene in cyclases]


def update(data: Dict[str, Any], record: SeqRecord) -> None:
    start = data["cluster"]["loci"].get("start_coord", 1) - 1
    end = data["cluster"]["loci"].get("end_coord", len(record.seq))
    genes = data["cluster"].get("genes", {})
    name_mapping = {}

    accession = data["cluster"]["loci"]["accession"]
    if "_" not in accession:
        return
    date = record.annotations["date"]

    extras = {}
    for extra in genes.get("extra_genes", []):
        assert extra["id"] not in extras
        exons = []
        if "location" not in extra:
            continue
        for exon in extra["location"]["exons"]:
            exons.append(FeatureLocation(exon["start"] - 1, exon["end"], extra["location"]["strand"]))
        if len(exons) > 1:
            location = CompoundLocation(exons)
        else:
            location = exons[0]
        extras[extra["id"]] = (location, extra)

    contained = []
    for feature in record.features:
        if feature.type == "CDS" and feature.location.start >= start and feature.location.end <= end:
            contained.append(feature)

    contained_by_name: Dict[str, SeqFeature] = {}
    duplicate_names = set()
    contained_locations = []
    for feature in contained:
        for qual in ["locus_tag", "gene", "protein_id"]:
            name = feature.qualifiers.get(qual, [""])[0]
            if not name:
                continue
            if contained_by_name.get(name, feature) != feature:
                duplicate_names.add(name)
                if "locus_tag" in feature.qualifiers and name != feature.qualifiers.get("locus_tag", [""])[0]:
                    continue
                raise ValueError(f"{data['cluster']['mibig_accession']}: multiple existing genes share a name: {name}, {feature.location}")
            contained_by_name[name] = feature
        contained_locations.append(feature.location)

    existing_annotations = {}
    for annotation in genes.get("annotations", []):
        existing_annotations[annotation["id"]] = annotation


    # finally, get to adding the information from the genbank
    for gene in contained:
        feature = None
        names = []
        already_present = False
        for qual in ["locus_tag", "protein_id", "gene"]:
            name = gene.qualifiers.get(qual, [""])[0]
            if name:
                names.append(name)
                if name in extras:
                    already_present = True

        for name in names:
            name_mapping[name] = names[0]
            if name in existing_annotations:
                existing_annotations[name]["id"] = names[0]

        if already_present:
            continue

        assert names, gene
        if "extra_genes" not in genes:
            genes["extra_genes"] = []
        genes["extra_genes"].append({
            "id": names[0],
            "location": {
                "exons": [{"end": part.end, "start": part.start + 1} for part in gene.location.parts],
                "strand": gene.location.strand,
            },
            "translation": gene.qualifiers["translation"][0],
        })
        annotation = {
            "id": names[0],
        }
        gene_name = gene.qualifiers.get("gene", [""])[0]
        if gene_name and annotation["id"] != gene_name and gene_name not in duplicate_names:
            annotation["name"] = gene_name
        if "product" in gene.qualifiers:
            assert '"' not in gene.qualifiers["product"]
            annotation["product"] = gene.qualifiers["product"][0]

        for name in names:
            if name not in existing_annotations:
                continue
            if "product" not in existing_annotations[name] and annotation.get("product"):
                existing_annotations[name]["product"] = annotation["product"]
            break
        else:
            if "annotations" not in genes:
                if "genes" not in data["cluster"]:
                    data["cluster"]["genes"] = genes
                genes["annotations"] = []
            genes["annotations"].append(annotation)

    update_gene_naming(data, name_mapping)

    comments = []
    if contained:
        comments.append(f"Gene details extracted from RefSeq accession {accession} dated {date}")
    if accession.startswith("NZ_"):
        data["cluster"]["loci"]["accession"] = data["cluster"]["loci"]["accession"][3:]
        comments.append("Changed from RefSeq to GenBank accession for stability")
    else:
        # try a manual mapping first
        new = FULL_ACCESSIONS.get(data["cluster"]["mibig_accession"])
        if new is None:
            # get the base accession sometimes without version from the comments
            new = get_refseq_base(record.annotations["comment"])
        assert new
        data["cluster"]["loci"]["accession"] = new
        
        comments.append("Changed from RefSeq to GenBank accession for stability")
    update_changelog(data, comments)


def update_changelog(data: Dict[str, Any], comments: List[str], contributor: str = "AAAAAAAAAAAAAAAAAAAAAAAA", release: str = "2.1") -> None:
    previous = data["changelog"][-1]

    if previous["version"] != release:
        new_entry = {
            "comments": comments,
            "contributors": [contributor],
            "version": release,
        }
        data["changelog"].append(new_entry)
        return
    for comment in comments:
        if comment not in previous["comments"]:
            previous["comments"].append(comment)
    if contributor not in previous["contributors"]:
        previous["contributors"].append(contributor)


def main(entry_file: str, gbk: str) -> int:
    with open(entry_file) as handle:
        data = json.load(handle)
    with open(gbk) as handle:
        records = list(SeqIO.parse(gbk, "genbank"))
    if len(records) != 1:
        print(f"too many records in genbank: {gbk}", file=sys.stderr)
        return 1
    try:
        empty, date = EMPTIES.get(data["cluster"]["mibig_accession"], (None, None))
        if empty:
            original = get_refseq_base(records[0].annotations["comment"])
            data["cluster"]["accession"] = original
            update_changelog(data, ["Changed from Refseq accession to equivalent GenBank accession"])
        else:
            update(data, records[0])
        with open(entry_file, "w") as handle:
            json.dump(data, handle, indent=4, separators=(',', ': '), sort_keys=True, ensure_ascii=False)
    except Exception:
        print(entry_file)
        traceback.print_exception(*sys.exc_info())
        return 1
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage:", os.path.basename(sys.argv[0]), "JSON GBK", file=sys.stderr)
        sys.exit(1)
    sys.exit(main(sys.argv[1], sys.argv[2]))
