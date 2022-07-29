#!/bin/bash

DATADIR=/gymreklab-tscc/helia/ensembl

bcftools view -h ${DATADIR}/ensemble_out/merged_chr1_sorted_ver2.vcf.gz > codis.vcf

for i in {1..22}
do
    bcftools view -H -R codis_loci.bed --no-update \
	${DATADIR}/ensemble_out/merged_chr"$i"_sorted_ver2.vcf.gz >> codis.vcf
done

bgzip -f codis.vcf
tabix -f -p vcf codis.vcf.gz

for pop in AFR AMR EUR SAS EAS
do
    statSTR \
	--vcf codis.vcf.gz \
	--vcftype hipstr \
	--afreq \
	--out freqs/"$pop" \
	--samples ../metadata/${pop}_names.txt
done
