#!/bin/bash
#PBS -q hotel 
#PBS -N mendelian
#PBS -l nodes=1:ppn=1
#PBS -l walltime=20:00:00
#PBS -o /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/mend/log
#PBS -e /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/mend/log
#PBS -V
#PBS -M hziaeija@ucsd.edu
#PBS -A gymreklab-group


cd /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/mend

source /projects/ps-gymreklab/helia/ensembl/venv_3.9/bin/activate

./MI_gt_analysis.py /projects/ps-gymreklab/helia/TR_1000G/1000G.ped chr${chrom} /projects/ps-gymreklab/helia/ensembl/ensemble_out/ensemble_chr${chrom}_filtered.vcf.gz > /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/mend/result/mend_gt_chr${chrom}.tab


#./MI_gt_analysis.py /projects/ps-gymreklab/helia/TR_1000G/1000G.ped chr${chrom} /projects/ps-gymreklab/helia/ensembl/1000G_calls/adVNTR/chrs/all_pop_filtered_chr${chrom}.vcf.gz > /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/mend/result/mend_gt_chr${chrom}_advntr.tab

