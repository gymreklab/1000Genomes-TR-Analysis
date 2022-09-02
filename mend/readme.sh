#!/bin/bash

chrom=21
./MI-analysis.py /gymreklab-tscc/helia/TR_1000G/1000G.ped chr${chrom} /gymreklab-tscc/helia/ensembl/ensemble_out/merged_chr${chrom}_sorted_ver2.vcf.gz > results/mend_chr${chrom}.tab
