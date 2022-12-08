#!/bin/bash

# Set up sample lists
KGFILE=/gymreklab-tscc/helia/TR_1000G/1000G.ped
cat ${KGFILE} | cut -d' ' -f 2 | grep -v Sample > all_samples.txt
for pop in EUR EAS SAS AMR AFR
do
    cat ${KGFILE} | grep -w ${pop} | cut -d' ' -f 2 > ${pop}_samples.txt
done
cat /gymreklab-tscc/helia/H3Africa/names/H3A_Baylor_sample_country.txt | \
    cut -f 1 > h3africa_samples.txt

# Compute het, n_alleles for all and for each superpopulation
DATADIR=/gymreklab-tscc/helia/ensembl/ensemble_out
for chrom in $(seq 1 22)
do
    echo $chrom
    statSTR \
	--vcf ${DATADIR}/ensemble_chr${chrom}_filtered.vcf.gz \
	--vcftype hipstr \
	--het --nalleles --nalleles-thresh 0.05 \
	--samples all_samples.txt,EUR_samples.txt,EAS_samples.txt,SAS_samples.txt,AMR_samples.txt,AFR_samples.txt,h3africa_samples.txt \
	--sample-prefixes ALL,EUR,EAS,SAS,AMR,AFR,H3A \
	--out chr${chrom}.stats
    bcftools query -f"%CHROM\t%POS\t%INFO/PERIOD\n" \
	${DATADIR}/ensemble_chr${chrom}_filtered.vcf.gz \
	> chr${chrom}.period.tab
done
