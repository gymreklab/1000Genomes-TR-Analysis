#!/bin/bash
#PBS -q hotel 
#PBS -N mendelian
#PBS -l nodes=1:ppn=1
#PBS -l walltime=20:00:00
#PBS -o /projects/ps-gymreklab/helia/ensembl/job_output
#PBS -e /projects/ps-gymreklab/helia/ensembl/job_output
#PBS -V
#PBS -M hziaeija@ucsd.edu
#PBS -A gymreklab-group


cd /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/mend

source /projects/ps-gymreklab/helia/ensembl/venv_3.9/bin/activate

./MI_gt_analysis.py /projects/ps-gymreklab/helia/TR_1000G/1000G.ped chr${chrom} /projects/ps-gymreklab/helia/ensembl/ensemble_out/merged_chr${chrom}_sorted_ver2.vcf.gz > /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/mend/result/mend_gt_chr${chrom}.tab
