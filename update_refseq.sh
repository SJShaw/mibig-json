#!/bin/bash

for bgc in `git grep '"accession"' | grep BGC | grep _ | grep -v retired | cut -d. -f1`; do
    if [ ! -f `basename $bgc`.gbk ] ; then
        wget https://mibig.secondarymetabolites.org/repository/`basename $bgc`/`basename $bgc`.gbk 2>/dev/null
        if [ "$?" != "0" ] ; then
            echo "no previous file fetchable for $bgc.json"
            continue
        fi 
    fi
    test -z "`grep -e 'Gene details extracted from RefSeq' -e 'equivalent GenBank accession' $bgc.json`"
    if [ "$?" == "0" ] ; then
        ./extract_refseq_info.py $bgc.json `basename $bgc`.gbk && \
        git commit $bgc.json -m "`basename $bgc`: add RefSeq gene details and change to GenBank accession"
        if [ "$?" != "0" ] ; then
            echo $bgc
            break
        fi
    fi
done
