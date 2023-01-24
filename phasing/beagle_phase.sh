#!/bin/bash

set -e

CHR=$1 # Target chromosome
REGIONS=$2 # File address containing TR coordinates one per line
SHORT_WINDOW_SIZE=$3

# Unphased TRs
VCF=/projects/ps-gymreklab/helia/ensembl/ensemble_out/ensemble_chr"$CHR"_filtered.vcf.gz

#Phased SNPs
SNPS=/projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/phasing/debug/1kGP_high_coverage_Illumina.chr"$CHR".filtered.SNV_INDEL_SV_phased_panel.vcf.gz

SAMPLES=/projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/phasing/samples.txt

echo "." > str_id

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
  SHORT_WINDOW="${chr}:$(( ${start} - ${SHORT_WINDOW_SIZE} ))-$(( ${end} + ${SHORT_WINDOW_SIZE} ))"
  echo $SHORT_WINDOW

  # Get unphased STRs
  bcftools view -S ${SAMPLES} -h $VCF --no-update > ensemble.${ID}.vcf
  bcftools view -S ${SAMPLES} -H $VCF --no-update -r $REGION | grep $start >> ensemble.${ID}.vcf
  bgzip -f ensemble.${ID}.vcf
  bcftools index -f ensemble.${ID}.vcf.gz

  # Subset samples
  bcftools view -S ${SAMPLES} -r $WINDOW $SNPS --no-update --output-type z --output-file snps.${ID}.vcf.gz
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
  java -Xmx5g -jar /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/phasing/beagle.r1399.jar \
    gt=snpstr.${ID}.fixmissing.vcf.gz \
    out="${ID}".phased usephase=true \
    nthreads=10 ped=pedigree.fam &&

  tabix -p vcf -f ${ID}.phased.vcf.gz

   # Select SNPs that close to STRs to fix phasing swap
    bcftools view -S ${SAMPLES} -r ${SHORT_WINDOW} ${ID}.phased.vcf.gz --no-update --output-type z --output-file ${ID}_truncated.phased.vcf.gz
    tabix -p vcf -f ${ID}_truncated.phased.vcf.gz


  # Make sure phase is consistent with original SNP file
    ./fix_switch_error.py \
        --phased-vcf "${ID}"_truncated.phased.vcf.gz \
        --ref-vcf $SNPS \
        --switch-threshold 0.5 --min-maf 0.1 --check-snps 500 \
        --new-vcf ${ID}_phased_fixed.vcf.gz &&
        tabix -p vcf ${ID}_phased_fixed.vcf.gz

  #extract STR only
    bcftools view -S ${SAMPLES} -i ID=@str_id ${ID}_phased_fixed.vcf.gz --no-update --output-type z --output-file ${ID}_str_phased_fixed.vcf.gz
    tabix -p vcf ${ID}_str_phased_fixed.vcf.gz
#  #remove unnecessary files
  rm ensemble.${ID}.vcf.gz* snps.${ID}.vcf.gz* snpstr.${ID}.vcf.gz* snpstr.${ID}.fixmissing.vcf.gz* ${ID}_truncated.phased.vcf.gz*

done <$REGIONS
