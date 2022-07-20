#!/bin/bash

method=$1
files=$2

for vcf in $files
do
    bcftools view -H -S samples.txt -R capillary_loci.bed $vcf --no-update --force-samples 
done > "$method"_cap_calls.txt
