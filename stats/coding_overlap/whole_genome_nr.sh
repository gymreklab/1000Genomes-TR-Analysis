#!/bin/bash

#PBS -N sample_nr
#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00
#PBS -o /projects/ps-gymreklab/helia/ensembl/experiments/stats/log
#PBS -e /projects/ps-gymreklab/helia/ensembl/experiments/stats/log
#PBS -A gymreklab-group
#PBS -q hotel
#PBS -V

# Getting alt GT for each sample
cd /projects/ps-gymreklab/helia/ensembl/experiments/stats

bcftools query -f '[%CHROM\t%POS\t%SAMPLE\t%GT\n]' -i 'GT="alt"' \
 /projects/ps-gymreklab/helia/ensembl/filtered_calls/chr"$chr"_MI_filtered.vcf.gz \
 > info/gt_chr"$chr".txt

