#!/bin/bash

#PBS -N motif-info
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -o /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/size-stats/log
#PBS -e /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/size-stats/log
#PBS -A gymreklab-group
#PBS -q hotel
#PBS -V

DATADIR=/projects/ps-gymreklab/helia/ensembl/ensemble_out
cd /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/het-stats

for chrom in $(seq 1 22); do bcftools query -f"%CHROM\t%POS\t%INFO/PERIOD\t%INFO/RU\n" ${DATADIR}/ensemble_chr${chrom}_filtered.vcf.gz; done > motif_info.tab
