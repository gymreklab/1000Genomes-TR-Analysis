#!/bin/bash


#PBS -N sample_nr
#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00
#PBS -o /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/log
#PBS -e /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/log
#PBS -A gymreklab-group
#PBS -q hotel
#PBS -V

# Getting GT for each sample

bcftools query -R /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/intersect.txt -f '[%CHROM\t%POS\t%SAMPLE\t%GT\n]' /projects/ps-gymreklab/helia/ensembl/ensemble_out/merged_chr"$chr"_sorted_ver2.vcf.gz > /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/info/gene_gt_chr"$chr".txt

