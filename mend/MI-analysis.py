#!/usr/bin/env python3
"""
Mendelian inheritance analysis for a single chrom

Usage: ./MI-anaylsis.py <pedigree> <chrom> <vcf>

Example:
./MI-analysis.py /gymreklab-tscc/helia/TR_1000G/1000G.ped chr21 /gymreklab-tscc/helia/ensembl/ensemble_out/merged_chr21_sorted_ver2.vcf.gz
"""

import pandas as pd
from cyvcf2 import VCF
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
from collections import defaultdict
import sys
import csv

try:
    pedfile = sys.argv[1]
    chrom = sys.argv[2]
    vcfpath = sys.argv[3]
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

pedigree = pd.read_csv(pedfile, delim_whitespace=True)
trios = pedigree[(pedigree['FatherID'] != "0") & (pedigree['MotherID'] != "0")]
child_in_trios = set(trios['SampleID'].to_list())
mother_in_trios = set(trios['MotherID'].to_list())
father_in_trios = set(trios['FatherID'].to_list())
all_ids = (child_in_trios.union(mother_in_trios)).union(father_in_trios)

IDs = trios[['SampleID', 'MotherID', 'FatherID']]
trio_IDs = []
for index, row in IDs.iterrows():
    trio_IDs.append((row['SampleID'], row['MotherID'], row['FatherID']))

sys.stderr.write("Found %s trios...\n"%(len(trio_IDs)))

def CheckMI(sample_GT, mother_GT, father_GT):
    if sample_GT[0] in mother_GT and sample_GT[1] in father_GT:
        return True
    if sample_GT[1] in mother_GT and sample_GT[0] in father_GT:
        return True
    return False

vcf = VCF(vcfpath, samples = list(all_ids))
samples = vcf.samples

for variant in vcf:
    n_families = 0  # Number of families with full call
    if len(variant.ALT) == 0: # If there is no alt allele here
        continue
    for family in trio_IDs:
        if "HG02567" in family: continue
        sample_index = samples.index(family[0])
        mother_index = samples.index(family[1])
        father_index = samples.index(family[2])
        fam_indices = [sample_index, mother_index, father_index]
        sample_GT = variant.genotypes[sample_index]
        mother_GT = variant.genotypes[mother_index]
        father_GT = variant.genotypes[father_index]
        if sample_GT[0] == -1 or mother_GT[0] == -1 or father_GT[0] == -1: # No call
            continue
        n_families += 1
        MI_val = CheckMI(sample_GT, mother_GT, father_GT)
        min_score = np.min([variant.format('SCOREGT')[ind] for ind in fam_indices])
        items = [variant.CHROM, variant.POS, variant.INFO["RU"], family[0], variant.INFO["METHODS"], MI_val, min_score]
        sys.stdout.write("\t".join([str(item) for item in items])+"\n")
        """
        score_GT_dict[int(variant.format('SCOREGT')[sample_index][0]*10)][0] += MI_val ## Follows MI
        score_GT_dict[int(variant.format('SCOREGT')[sample_index][0]*10)][1] += 1 ## Total number
        loci_dict[variant.POS] += MI_val
        """
