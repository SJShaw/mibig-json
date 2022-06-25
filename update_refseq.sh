#!/bin/bash

for bgc in `git grep '"accession"' | grep BGC | grep _ | grep -v retired | cut -d. -f1`; do
    if [ ! -f `basename $bgc`.gbk ] ; then
        wget https://mibig.secondarymetabolites.org/repository/$bgc/$bgc.gbk 2>/dev/null
        if [ "$?" != "0" ] ; then
            echo "no previous file fetchable for $bgc.json"
            break
        fi 
    fi
    test -z "`grep 'Gene details extracted from RefSeq' $bgc.json`"
    if [ "$?" == "0" ] ; then
        ./extract_refseq_info.py $bgc.json `basename $bgc`.gbk && \
        git commit $bgc.json -m "`basename $bgc`: add RefSeq gene details and change to GenBank accession"
        if [ "$?" != "0" ] ; then
            break
        fi
    fi
done
