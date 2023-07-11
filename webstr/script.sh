#!/bin/bash

#PBS -N webstr
#PBS -l nodes=1:ppn=1
#PBS -l walltime=3:00:00
#PBS -o /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/webstr/log
#PBS -e /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/webstr/log
#PBS -A gymreklab-group
#PBS -q hotel
#PBS -V

source /home/hziaeija/miniconda3/bin/activate
conda init bash
conda activate ensembl-vep

webstr_addr="/projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/webstr"

vep -i /projects/ps-gymreklab/helia/ensembl/ensemble_out/ensemble_chr"$chr"_filtered.vcf.gz \
    -o "$webstr_addr"/genome_annotations/gene_annotation_chr${chr}.txt \
    --force_overwrite -gtf /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/webstr/genome_annotations/gencode.v22.annotation.gtf.gz \
    --fasta /projects/ps-gymreklab/helia/ensembl/hg38.analysisSet.fa


bcftools query -f '%CHROM\t%POS\t%END\t%RU\n' /projects/ps-gymreklab/helia/ensembl/ensemble_out/ensemble_chr"$chr"_filtered.vcf.gz > "$webstr_addr"/tables/repeat_tables/ensemble_"$chr".txt


python3 "$webstr_addr"/create_tables.py $chr

