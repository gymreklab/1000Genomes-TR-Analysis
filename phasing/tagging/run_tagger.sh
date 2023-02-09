#!/bin/bash


REGIONS=$1
pop=$2
chr=$3

SNP_STRs=/projects/ps-gymreklab/yal084_storage/ensemble_phasing/chr"$3"_phased_combined_sorted.vcf.gz
SNPs=/projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/phasing/snps/1kGP_high_coverage_Illumina.chr"$3".filtered.SNV_INDEL_SV_phased_panel.vcf.gz

while read p; do
  chr=$(echo $p | cut -d " " -f 1)
  start=$(echo $p | cut -d " " -f 2)
  end=$(echo $p | cut -d " " -f 3)

  window_start=$(($start - 50000))
  window_end=$(($end + 50000))
  REGION="$chr:$start-$end"
  WINDOW="$chr:$window_start-$window_end"


  bcftools view -H -r $REGION -S /projects/ps-gymreklab/helia/ensembl/experiments/allele_freq/names/"$pop"_names.txt --no-update $SNP_STRs | grep $start > "$REGION"_"$pop".txt &&

  bcftools view -H -r $WINDOW -S /projects/ps-gymreklab/helia/ensembl/experiments/allele_freq/names/"$pop"_names.txt --no-update $SNPs > "$WINDOW"_"$pop".txt &&

  python3 tagged_strs.py "$REGION"_"$pop".txt "$WINDOW"_"$pop".txt &&

  rm "$REGION"_"$pop".txt "$WINDOW"_"$pop".txt

done < $REGIONS
