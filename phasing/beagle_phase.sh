#!/bin/bash

set -e

CHR=$1 # Target chromosome
REGIONS=$2 # File address containing TR coordinates one per line

# Unphased TRs
VCF=/projects/ps-gymreklab/helia/ensembl/phasing_st/filtered/chr"$CHR"_MI_filtered.vcf.gz

#Phased SNPs
SNPS=/projects/ps-gymreklab/helia/ensembl/biallelic/snps/CCDG_14151_B01_GRM_WGS_2020-08-05_chr"$CHR".filtered.shapeit2-duohmm-phased.vcf.gz

while read p; do
  # Extracting region and window
  chr=$(echo $p | cut -d " " -f 1)
  start=$(echo $p | cut -d " " -f 2)
  end=$(echo $p | cut -d " " -f 3)
  window_start=$(($start - 50000))
  window_end=$(($end + 50000))
  REGION="$chr:$start-$end"
  WINDOW="$chr:$window_start-$window_end"
  ID="$chr-$start-$end"
  echo $ID

  # Get unphased STRs
  bcftools view -S samples.txt -r $REGION $VCF --no-update --output-type z --output-file ensemble.${ID}.vcf.gz
  bcftools index -f ensemble.${ID}.vcf.gz

  # Subset samples
  bcftools view -S samples.txt -r $WINDOW $SNPS --no-update --output-type z --output-file snps.${ID}.vcf.gz
  bcftools index -f snps.${ID}.vcf.gz

  # Concatenate STRs+SNPs
  bcftools concat -a snps.${ID}.vcf.gz ensemble.${ID}.vcf.gz -O z -o snpstr.${ID}.vcf.gz &&
  bcftools index -f snpstr.${ID}.vcf.gz

  # Fix the missing marker
  bcftools view snpstr.${ID}.vcf.gz |
    sed 's/\.\:\./\.\/\.\:\./g' |
    sed 's/\.\/\.\/\./\.\/\./g' | bgzip -c > snpstr.${ID}.fixmissing.vcf.gz &&
    bcftools index -f -t snpstr.${ID}.fixmissing.vcf.gz

  # Phase with Beagle
  java -Xmx5g -jar /projects/ps-gymreklab/helia/ensembl/phasing_st/phased/beagle.r1399.jar \
    gt=snpstr.${ID}.fixmissing.vcf.gz \
    out="${ID}".phased usephase=true \
    nthreads=4 &&

  tabix -p vcf -f "${ID}".phased.vcf.gz

  # Make sure phase is consistent with original SNP file
  ./fix_switch_error.py \
    --phased-vcf "${ID}".phased.vcf.gz \
    --ref-vcf snps.${ID}.vcf.gz \
    --switch-threshold 0.5 --min-maf 0.1 --check-snps 100 \
    --new-vcf beagle_out/"$CHR"/${ID}_phased_fixed.vcf.gz &&

  tabix -p vcf beagle_out/"$CHR"/${ID}_phased_fixed.vcf.gz

  #remove unnecessary files
  rm ensemble.${ID}.vcf.gz* snps.${ID}.vcf.gz* snpstr.${ID}.vcf.gz* snpstr.${ID}.fixmissing.vcf.gz* "${ID}".phased.vcf.gz* "${ID}".phased.log

done <$REGIONS
