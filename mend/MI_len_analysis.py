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

def CheckMI(sample_alleles, mother_alleles, father_alleles):
    sample_alleles = [len(al) for al in sample_alleles.split("/")]
    mother_alleles = [len(al) for al in mother_alleles.split("/")]
    father_alleles = [len(al) for al in father_alleles.split("/")]
    if sample_alleles[0] in mother_alleles and sample_alleles[1] in father_alleles:
        return True
    if sample_alleles[1] in mother_alleles and sample_alleles[0] in father_alleles:
        return True
    return False

def IsRef(gt):
    return gt[0]==0 and gt[1]==0

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
        if IsRef(sample_GT) and IsRef(mother_GT) and IsRef(father_GT): # all homozygous ref
            continue
        #gbs=("%s,%s,%s"%(variant.format("GB")[sample_index],
        #                  variant.format("GB")[mother_index],
        #                  variant.format("GB")[father_index]))
        MI_val = CheckMI(variant.gt_bases[sample_index], variant.gt_bases[mother_index], variant.gt_bases[father_index])
        min_score_gt = np.min([variant.format('SCORE')[ind] for ind in fam_indices])
        items = [variant.CHROM, variant.POS, variant.INFO["RU"], family[0], \
                 MI_val, min_score_gt]
        sys.stdout.write("\t".join([str(item) for item in items])+"\n")

