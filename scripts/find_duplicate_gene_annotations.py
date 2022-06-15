#!/usr/bin/env python3

import glob
import json
import os
import sys


def run_all():
    for dir_name in ["data", "pending"]:
        for file in glob.glob(os.path.join(dir_name, "*.json")):
            run_single(file, prefix=dir_name + os.sep)


def run_single(file: str, prefix: str = ""):
    with open(file) as handle:
        try:
            data = json.load(handle)
        except json.JSONDecodeError as err:
            print(f"{prefix}{os.path.split(file)[-1]}: {err}")
            return
    try:
        ids = []
        names = []
        try:
            for gene in data["cluster"].get("genes", {}).get("annotations", {}):
                ids.append(gene["id"])
                name = gene.get("name")
                if name:
                    names.append(name)
        except KeyError:
            print(os.path.split(file)[-1])
            raise
        id_mismatch = len(set(ids)) != len(ids)
        name_mismatch = len(set(names)) != len(names)
        if id_mismatch or name_mismatch:
            print(f"{prefix}{os.path.split(file)[-1]}: duplicate gene annotations for {'ids' if id_mismatch else ''}{' and ' if id_mismatch and name_mismatch else ''}{'names' if name_mismatch else ''}")
    except ValidationError as err:
        lines = str(err).splitlines()
        print(f"{prefix}{os.path.split(file)[-1]}: {lines[-2]} {lines[0]}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        run_all()
    else:
        for filename in sys.argv[1:]:
            run_single(filename)
