#!/bin/bash
#PBS -q hotel 
#PBS -N mendelian_gangstr
#PBS -l nodes=1:ppn=1
#PBS -l walltime=40:00:00
#PBS -o /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/mend/log
#PBS -e /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/mend/log
#PBS -V
#PBS -M hziaeija@ucsd.edu
#PBS -A gymreklab-group


cd /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/mend

source /projects/ps-gymreklab/helia/ensembl/venv_3.9/bin/activate

./MI_gt_analysis.py /projects/ps-gymreklab/helia/TR_1000G/1000G.ped chr${chrom} /projects/ps-gymreklab/helia/ensembl/ensemble_out/ensemble_chr${chrom}_filtered.vcf.gz > /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/mend/result/mend_gt_chr${chrom}_extra_info_r2.tab


#./MI_gt_analysis.py /projects/ps-gymreklab/helia/TR_1000G/1000G.ped chr${chrom} /home/saraj/vntr/ensembleTR/filtered_calls/filtered_advntr_calls_str_like_filter/all_pop_filtered_plus_str_like_filter_chr"$chrom".vcf.gz > /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/mend/result/mend_gt_chr${chrom}_advntr.tab


#./MI_gt_analysis.py /projects/ps-gymreklab/helia/TR_1000G/1000G.ped chr${chrom} /projects/ps-gymreklab/helia/ensembl/1000G_calls/gangstr/chrs/gangstr_all_merged_chr${chrom}_corrected.vcf.gz > /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/mend/result/mend_gt_chr${chrom}_gangstr.tab

#./MI_len_analysis_hipstr.py /projects/ps-gymreklab/helia/TR_1000G/1000G.ped chr${chrom} /projects/ps-gymreklab/helia/ensembl/1000G_calls/hipstr/chrs/hipstr_corrected_chr${chrom}.vcf.gz > /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/mend/result/mend_len_chr${chrom}_hipstr.tab

#./MI_len_analysis.py /projects/ps-gymreklab/helia/TR_1000G/1000G.ped chr${chrom} /projects/ps-gymreklab/helia/ensembl/ensemble_out/before_revision_files/ensemble_chr${chrom}_filtered.vcf.gz > /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/mend/result/mend_len_chr${chrom}.tab
