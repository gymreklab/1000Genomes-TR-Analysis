#!/bin/bash


REGIONS=$1

while read p; do
  chr=$(echo $p | cut -d " " -f 1)
  start=$(echo $p | cut -d " " -f 2)
  end=$(echo $p | cut -d " " -f 3)

  window_start=$(($start - 50000))
  window_end=$(($end + 50000))
  REGION="$chr:$start-$end"
  WINDOW="$chr:$window_start-$window_end"

  SNPs=/projects/ps-gymreklab/helia/ensembl/biallelic/snps/CCDG_14151_B01_GRM_WGS_2020-08-05_"$chr".filtered.shapeit2-duohmm-phased.vcf.gz
  SNP_STRs=/projects/ps-gymreklab/helia/ensembl/phasing_st/phased/beagle_out/21/merged/"$chr".phased.vcf.gz

  bcftools view  -H -r $REGION -S samples.txt --no-update $SNP_STRs > "$REGION".txt &&

  bcftools view -S samples.txt -r $WINDOW $SNPs --no-update -H > "$WINDOW".txt &&

  python3 tagged_strs.py "$REGION".txt "$WINDOW".txt &&

  rm "$REGION".txt "$WINDOW".txt

done < $REGIONS
