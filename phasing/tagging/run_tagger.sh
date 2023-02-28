#!/bin/bash


REGIONS=$1
pop=$2
chr=$3

STRs=/projects/ps-gymreklab/yal084_storage/ensemble_phasing/annotation_fixed_no_snp/chr"$i"_SNP_merged_final.vcf.gz
SNPs=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr"$chr".filtered.SNV_INDEL_SV_phased_panel.vcf.gz

while read p; do
  chr=$(echo $p | cut -d " " -f 1)
  start=$(echo $p | cut -d " " -f 2)
  end=$(echo $p | cut -d " " -f 3)

  window_start=$(($start - 50000))
  window_end=$(($end + 50000))
  REGION="$chr:$start-$end"
  WINDOW="$chr:$window_start-$window_end"

  bcftools view -H -r $REGION -S /projects/ps-gymreklab/helia/ensembl/experiments/allele_freq/names/"$pop"_names.txt --no-update $STRs | grep $start > "$REGION"_"$pop".txt

  bcftools view -H -r $WINDOW -S /projects/ps-gymreklab/helia/ensembl/experiments/allele_freq/names/"$pop"_names.txt --no-update $SNPs > "$WINDOW"_"$pop".txt

  python3 tagged_strs.py "$REGION"_"$pop".txt "$WINDOW"_"$pop".txt

#  time python3 tagged_strs_alt.py STR.vcf.gz SNPS.vcf

  rm "$REGION"_"$pop".txt "$WINDOW"_"$pop".txt

done < $REGIONS
