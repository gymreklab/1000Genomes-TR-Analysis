#!/bin/bash

#PBS -N dbsnp_fix
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -o /projects/ps-gymreklab/helia/ensembl/experiments/dbsnp/log
#PBS -e /projects/ps-gymreklab/helia/ensembl/experiments/dbsnp/log
#PBS -A gymreklab-group
#PBS -q hotel
#PBS -V

declare -A chrs=( ["1"]="NC_000001.11" ["2"]="NC_000002.12" ["3"]="NC_000003.12" ["4"]="NC_000004.12" ["5"]="NC_000005.10" ["6"]="NC_000006.12" ["7"]="NC_000007.14" ["8"]="NC_000008.11" ["9"]="NC_000009.12" ["10"]="NC_000010.11" ["11"]="NC_000011.10" ["12"]="NC_000012.12" ["13"]="NC_000013.11" ["14"]="NC_000014.9" ["15"]="NC_000015.10" ["16"]="NC_000016.10" ["17"]="NC_000017.11" ["18"]="NC_000018.10" ["19"]="NC_000019.10" ["20"]="NC_000020.11" ["21"]="NC_000021.9" ["22"]="NC_000022.11")


echo $chr
echo ${chrs["$chr"]}
cd /projects/ps-gymreklab/helia/ensembl/experiments/dbsnp

bcftools view -r ${chrs["$chr"]} GCF_000001405.39.gz --no-update -O z -o snps/dbsnp_hg38_"$chr".vcf.gz &&

tabix -p vcf snps/dbsnp_hg38_"$chr".vcf.gz &&

bcftools norm -m-any snps/dbsnp_hg38_"$chr".vcf.gz -O z -o snps/dbsnp_hg38_"$chr"_biallelic.vcf.gz &&

tabix -p vcf snps/dbsnp_hg38_"$chr"_biallelic.vcf.gz


