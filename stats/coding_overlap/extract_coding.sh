#!/bin/bash


#PBS -N extract_coding
#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00
#PBS -o /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/log
#PBS -e /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/log
#PBS -A gymreklab-group
#PBS -q hotel
#PBS -V


source /projects/ps-gymreklab/helia/ensembl/venv_3.9/bin/activate
FILE_DIR=/projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/stats/coding_overlap
OUT_DIR=/projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/heterozygosity
cd $FILE_DIR

bcftools view -h /projects/ps-gymreklab/helia/ensembl/filtered_calls/chr"$chr"_MI_filtered.vcf.gz > \
$OUT_DIR/coding_calls_"$chr".vcf &&

grep -w chr"$chr" $FILE_DIR/TR_intersect.txt | awk '{print $2}' > $FILE_DIR/pos_chr"$chr".txt &&

bcftools view -H /projects/ps-gymreklab/helia/ensembl/filtered_calls/chr"$chr"_MI_filtered.vcf.gz | grep -F -f $FILE_DIR/pos_chr"$chr".txt >> $OUT_DIR/coding_calls_"$chr".vcf

bgzip -f $OUT_DIR/coding_calls_"$chr".vcf
tabix -p vcf -f $OUT_DIR/coding_calls_"$chr".vcf.gz

rm pos_chr"$chr".txt
