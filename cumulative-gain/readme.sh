#!/bin/bash

# Sample info 1000G: /gymreklab-tscc/helia/TR_1000G/1000G.ped
# Sample info H3Africa: /gymreklab-tscc/helia/H3Africa/names/H3A_Baylor_sample_country.txt

# VCFs: 
# /gymreklab-tscc/helia/ensembl/ensemble_out/ensemble_chr*_filtered.vcf.gz

# Compute only on chr1, should be representative
# x-axis: add one sample at a time. y-axis: number of alleles
# Exclude homopolymers. break into: singleton, doubleton, polymorphic
# Order samples: EUR, EAS, AMR, SAS, AFR, H3Africa
./get_cumulative_counts.py > cumulative_allele_counts.tab

# Plot with notebook
