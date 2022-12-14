#!/bin/bash

CHAINFILE=/storage/resources/dbase/human/hg38/hg19ToHg38.over.chain.gz

for f in $(ls raw/master/*.tab)
do 
    tissue=$(basename $f _master.tab)
    echo $tissue
    infile=raw/master/${tissue}_master.tab
    cat $infile | grep -v mashr | awk -F "\t" '{print $1 "\t" $3 "\t" $3+1 "\t" $1":"$3}' \
	> liftover/${tissue}_hg19.bed
    liftOver liftover/${tissue}_hg19.bed ${CHAINFILE} liftover/${tissue}_hg38.bed liftover/${tissue}_unMapped
done
