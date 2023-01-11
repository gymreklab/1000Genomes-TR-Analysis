#!/bin/bash

PREFIX=$1
DATADIR=/projects/ps-gymreklab/helia/ensembl
#DATADIR=/projects/ps-gymreklab/helia/ensembl

bcftools view -h ${DATADIR}/ensemble_out/ensemble_chr1_filtered.vcf.gz  > ${PREFIX}.vcf

for i in {1..22}
do
    bcftools view -H -R problematic_loci.txt --no-update \
        ${DATADIR}/ensemble_out/ensemble_chr"$i"_filtered.vcf.gz  >> ${PREFIX}.vcf
done

bgzip -f ${PREFIX}.vcf
tabix -f -p vcf ${PREFIX}.vcf.gz


statSTR \
	--vcf ${PREFIX}.vcf.gz \
        --vcftype hipstr \
        --afreq \
        --out "$PREFIX"
