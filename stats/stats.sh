#!/bin/bash

#PBS -N afreq
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:30:00
#PBS -o /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/log
#PBS -e /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/log
#PBS -A gymreklab-group
#PBS -q hotel
#PBS -V

# Getting statistics about genes and whole genome

cd /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions

bcftools query -R intersect.txt -f '%CHROM\t%POS\t%PERIOD\t%ALT\n' /projects/ps-gymreklab/helia/ensembl/ensemble_out/merged_chr"$chr"_sorted_ver2.vcf.gz > info/gene_info_chr"$chr".txt

bcftools query  -f '%CHROM\t%POS\t%PERIOD\t%ALT\n' /projects/ps-gymreklab/helia/ensembl/ensemble_out/merged_chr"$chr"_sorted_ver2.vcf.gz > info/all_info_chr"$chr".txt
