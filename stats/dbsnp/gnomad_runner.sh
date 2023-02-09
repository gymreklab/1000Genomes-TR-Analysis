#!/bin/bash

#PBS -N gnomad
#PBS -l nodes=1:ppn=5
#PBS -l walltime=30:00:00
#PBS -o /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/stats/dbsnp/log
#PBS -e /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/stats/dbsnp/log
#PBS -A gymreklab-group
#PBS -q hotel
#PBS -V



source /projects/ps-gymreklab/helia/ensembl/venv_3.9/bin/activate

cd /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/stats/dbsnp

echo $chr

python3 dbsnp_overlap.py /projects/ps-gymreklab/helia/ensembl/hg38.analysisSet.fa \
			 https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr"$chr".vcf.bgz \
       /projects/ps-gymreklab/helia/ensembl/ensemble_out/ensemble_chr"$chr"_filtered.vcf.gz \
       $chr gnomad
