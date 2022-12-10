#!/bin/bash

#PBS -N matrix
#PBS -l nodes=1:ppn=5
#PBS -l walltime=10:00:00
#PBS -o /projects/ps-gymreklab/helia/ensembl/experiments/pca/log
#PBS -e /projects/ps-gymreklab/helia/ensembl/experiments/pca/log
#PBS -V
#PBS -A gymreklab-group
#PBS -q hotel

source /projects/ps-gymreklab/helia/ensembl/venv_3.9/bin/activate


echo $chr

zcat /projects/ps-gymreklab/helia/ensembl/filtered_calls/African_calls/chr"$chr"_MI_filtered_AFR.vcf.gz | python3 /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/stats/pca/vcf_mat.py $chr

echo "Done"
