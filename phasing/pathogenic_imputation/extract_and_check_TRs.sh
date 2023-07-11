#!/bin/bash



for i in {1..22}
do	
	bcftools view -R blood_loci.txt --no-update /projects/ps-gymreklab/helia/ensembl/ensemble_out/ensemble_chr"$i"_filtered.vcf.gz -O z -o blood_vcf/LOO_chr"$i".vcf.gz
	tabix -f -p vcf blood_vcf/LOO_chr"$i".vcf.gz
	bcftools view -R /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/capillary/capillary_loci.bed /projects/ps-gymreklab/helia/ensembl/ensemble_out/ensemble_chr"$i"_filtered.vcf.gz -O z -o pathogenic_vcf/LOO_chr"$i".vcf.gz
	tabix -f -p vcf pathogenic_vcf/LOO_chr"$i".vcf.gz
done
