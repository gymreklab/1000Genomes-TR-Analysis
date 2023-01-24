#!/bin/bash

SNP_addr=$1
samples=/projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/phasing/samples.txt
REGIONS=$2
STR_addr=$3

while read p; do
  # Extracting region and window
  chr=$(echo $p | cut -d " " -f 1)
  start=$(echo $p | cut -d " " -f 2)
  end=$(echo $p | cut -d " " -f 3)
  window_start=$(($start - 50000))
  window_end=$(($end + 50000))
  REGION="$chr:$start-$end"
  WINDOW="$chr:$window_start-$window_end"

  bcftools view -S $samples -H -r $REGION --no-update $STR_addr > "$STR_addr".txt
  bcftools view -S $samples -H -r $WINDOW --no-update $SNP_addr > "$SNP_addr".txt


  python3 /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/phasing/tagging/tagged_strs.py "$STR_addr".txt "$SNP_addr".txt


  rm "$STR_addr".txt "$SNP_addr".txt

done < $REGIONS
