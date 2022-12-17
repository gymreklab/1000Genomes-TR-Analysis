#!/bin/bash


REGIONS=$1
pop=$2

while read p; do
  chr=$(echo $p | cut -d " " -f 1)
  start=$(echo $p | cut -d " " -f 2)
  end=$(echo $p | cut -d " " -f 3)

  window_start=$(($start - 50000))
  window_end=$(($end + 50000))
  REGION="$chr:$start-$end"
  WINDOW="$chr:$window_start-$window_end"

  SNP_STRs=/projects/ps-gymreklab/helia/ensembl/phasing_st/phased/beagle_out/21/merged/"$chr".phased.vcf.gz

  bcftools view -H -r $REGION -S /projects/ps-gymreklab/helia/ensembl/experiments/allele_freq/names/"$pop"_names.txt --no-update $SNP_STRs > "$REGION"_"$pop".txt &&

  bcftools view -H -r $WINDOW -S /projects/ps-gymreklab/helia/ensembl/experiments/allele_freq/names/"$pop"_names.txt --no-update $SNP_STRs > "$WINDOW"_"$pop".txt &&

  python3 tagged_strs.py "$REGION"_"$pop".txt "$WINDOW"_"$pop".txt &&

  rm "$REGION"_"$pop".txt "$WINDOW"_"$pop".txt

done < $REGIONS
