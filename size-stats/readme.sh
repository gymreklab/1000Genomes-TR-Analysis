#!/bin/bash

# For each sample, count how many alleles match ref, ref+1, ref-1, etc.
# Do separately for homopolymers and non-homopolymers
for chrom in $(seq 1 22)
do
    ./summarize_allele_sizes.py ${chrom} > asize_summary_chr${chrom}.tab
done

# Get list of STR, period, diff, AF
for chrom in $(seq 1 22)
do
    ./get_diff_af.py ${chrom}
done > diff_afreqs.tab
