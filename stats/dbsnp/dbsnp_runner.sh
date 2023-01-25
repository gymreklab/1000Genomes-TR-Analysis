#!/bin/bash

#PBS -N dbsnp
#PBS -l nodes=1:ppn=2
#PBS -l walltime=2:00:00
#PBS -o /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/stats/dbsnp/log
#PBS -e /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/stats/dbsnp/log
#PBS -A gymreklab-group
#PBS -q hotel



source /projects/ps-gymreklab/helia/ensembl/venv_3.9/bin/activate

cd /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/stats/dbsnp

python3 dbsnp_overlap.py /projects/ps-gymreklab/helia/ensembl/hg38.analysisSet.fa \
			 /projects/ps-gymreklab/helia/ensembl/experiments/dbsnp/snps/dbsnp_hg38_"$chr"_biallelic.vcf.gz \
       /projects/ps-gymreklab/helia/ensembl/ensemble_out/ensemble_chr"$chr"_filtered.vcf.gz \
       $chr
