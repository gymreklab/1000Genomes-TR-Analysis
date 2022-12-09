#!/bin/bash

PREFIX=$1
DATADIR=/projects/ps-gymreklab/helia/ensembl
#DATADIR=/projects/ps-gymreklab/helia/ensembl

bcftools view -h ${DATADIR}/filtered_calls/chr1_MI_filtered.vcf.gz  > ${PREFIX}.vcf

for i in {1..22}
do
    bcftools view -H -R ${PREFIX}_loci.bed --no-update \
	${DATADIR}/filtered_calls/chr"$i"_MI_filtered.vcf.gz  >> ${PREFIX}.vcf
done

bgzip -f ${PREFIX}.vcf
tabix -f -p vcf ${PREFIX}.vcf.gz

for pop in AFR AMR EUR SAS EAS H3Africa
do
    statSTR \
	--vcf ${PREFIX}.vcf.gz \
	--vcftype hipstr \
	--afreq \
	--out freqs/${PREFIX}-${pop} \
	--samples /projects/ps-gymreklab/helia/ensembl/experiments/allele_freq/names/${pop}_names.txt
done
