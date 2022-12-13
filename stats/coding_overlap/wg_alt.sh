#!/bin/bash


#PBS -N wg_polymorphism
#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00
#PBS -o /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/log
#PBS -e /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/log
#PBS -A gymreklab-group
#PBS -q hotel
#PBS -V


bcftools query -f '%CHROM\t%POS\t%PERIOD\t%END\t%ALT\n' /projects/ps-gymreklab/helia/ensembl/filtered_calls/chr"$chr"_MI_filtered.vcf.gz  > /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/info/stats_chr"$chr".txt

