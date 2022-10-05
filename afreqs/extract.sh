#!/bin/bash

PREFIX=$1
DATADIR=/gymreklab-tscc/helia/ensembl
#DATADIR=/projects/ps-gymreklab/helia/ensembl

bcftools view -h ${DATADIR}/ensemble_out/merged_chr1_sorted_ver2.vcf.gz > ${PREFIX}.vcf

for i in {1..22}
do
    bcftools view -H -S ../metadata/samples.txt -R ${PREFIX}_loci.bed --no-update \
	${DATADIR}/ensemble_out/merged_chr"$i"_sorted_ver2.vcf.gz >> ${PREFIX}.vcf
done

bgzip -f ${PREFIX}.vcf
tabix -f -p vcf ${PREFIX}.vcf.gz

for pop in AFR AMR EUR SAS EAS
do
    statSTR \
	--vcf ${PREFIX}.vcf.gz \
	--vcftype hipstr \
	--afreq \
	--out freqs/${PREFIX}-${pop} \
	--samples ../metadata/${pop}_names.txt
done
