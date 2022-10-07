#!/bin/bash


#PBS -N heterozygosity
#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00
#PBS -o /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/log
#PBS -e /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/log
#PBS -A gymreklab-group
#PBS -q hotel
#PBS -V


source /projects/ps-gymreklab/helia/ensembl/venv_3.9/bin/activate
cd /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/stats/coding_overlap
#bcftools view -h /projects/ps-gymreklab/helia/ensembl/ensemble_out/merged_chr"$chr"_sorted_ver2.vcf.gz > \
#/projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/heterozygosity/coding_calls_"$chr".vcf
#
#grep chr"$chr" TR_intersect.txt | awk '{print $2}' > pos_chr"$chr".txt
#
#bcftools view -H /projects/ps-gymreklab/helia/ensembl/ensemble_out/merged_chr"$chr"_sorted_ver2.vcf.gz | grep -F -f pos_chr"$chr".txt \
#>> /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/heterozygosity/coding_calls_"$chr".vcf
#
#
#bgzip -f /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/heterozygosity/coding_calls_"$chr".vcf
#tabix -p vcf -f /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/heterozygosity/coding_calls_"$chr".vcf.gz
#
#rm pos_chr"$chr".txt


time statSTR \
        --vcf /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/heterozygosity/coding_calls_"$chr".vcf.gz \
        --vcftype hipstr \
        --het \
        --out /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/heterozygosity/coding_het_"$pop_name"_"$chr" \
        --samples /projects/ps-gymreklab/helia/ensembl/experiments/allele_freq/names/"$pop_name"_names.txt

