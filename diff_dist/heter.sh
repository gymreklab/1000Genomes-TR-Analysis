#!/bin/bash


#PBS -N heterozygosity
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -o /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/log
#PBS -e /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/log
#PBS -A gymreklab-group
#PBS -q hotel
#PBS -V



source /projects/ps-gymreklab/helia/ensembl/venv_3.9/bin/activate

time statSTR \
        --vcf /projects/ps-gymreklab/helia/ensembl/ensemble_out/merged_chr"$chr"_sorted_ver2.vcf.gz \
        --vcftype hipstr \
        --afreq \
	--mean \
        --out /projects/ps-gymreklab/helia/ensembl/experiments/charact/diff_dist/heterozygosity/coding_het_"$pop_name"_"$chr" \
        --samples /projects/ps-gymreklab/helia/ensembl/experiments/allele_freq/names/"$pop_name"_names.txt
