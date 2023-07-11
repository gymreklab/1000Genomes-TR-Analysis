#!/usr/bin/env python3

import cyvcf2
import numpy as np
import pandas as pd
import sys

########## Load sample info and get sample order #############
kgsamples = pd.read_csv("/projects/ps-gymreklab/helia/TR_1000G/1000G.ped", delim_whitespace=True)
h3afrsamples = pd.read_csv("/projects/ps-gymreklab/helia/H3Africa/names/H3A_Baylor_sample_country.txt", \
                           sep="\t", names=["SampleID","Country"])
h3afrsamples["Superpopulation"] = "H3Africa"
h3afrsamples["Population"] = h3afrsamples["Country"]
samples = pd.concat([kgsamples[["SampleID","Population","Superpopulation"]], \
                     h3afrsamples[["SampleID","Population","Superpopulation"]]])

# Sort by correct order
def GetPopOrder(spop):
    spops = ["EUR","EAS","SAS","AMR","AFR","H3Africa"]
    return "%s_%s"%(spops.index(spop), spop)

samples["Superpopulation"] = samples["Superpopulation"].apply(GetPopOrder)
samples = samples.sort_values(["Superpopulation","Population"])
sample_list = list(samples["SampleID"])

############# Preprocess VCF ##########
# Get strid:allele -> [ordered list of samples]
# Update counts of sample->sing, doub, polymorphic

# Set up sample counts info
sample_counts = {}
for s in list(samples["SampleID"]): sample_counts[s] = {"single": 0, "double":0, "poly": 0}

maxloc = np.inf # for debugging
numloc = 0

# Go through VCF
reader = cyvcf2.VCF("/projects/ps-gymreklab/helia/ensembl/ensemble_out/ensemble_chr21_filtered.vcf.gz")
vcf_sample_list = reader.samples
sample_order = dict(zip(sample_list, [vcf_sample_list.index(s) for s in sample_list]))
sample_pop = dict(zip(list(samples["SampleID"]), list(samples["Population"])))

for v in reader:
    if v.FILTER is not None: continue
    if v.INFO["PERIOD"] == 1: continue
    acounts = {}
    sys.stderr.write("%s\n"%v.POS)
    for s in sample_list:
        gt = v.genotypes[sample_order[s]]
        for allele in gt[0:2]:
            if allele not in acounts:
                acounts[allele] = 1
            else:
                acounts[allele] = acounts[allele] + 1
        # counts are cumulative, based on including up to this sample
        for a in acounts:
            if acounts[a] == 1:
                sample_counts[s]["single"] = sample_counts[s]["single"] + 1
            elif acounts[a] == 2:
                sample_counts[s]["double"] = sample_counts[s]["double"] + 1
            else:
                sample_counts[s]["poly"] = sample_counts[s]["poly"] + 1
    numloc += 1
    if numloc == maxloc: break

##### Print out data ####
for s in sample_list:
    sys.stdout.write("\t".join([s, str(sample_counts[s]["single"]), str(sample_counts[s]["double"]), \
                                str(sample_counts[s]["poly"]), sample_pop[s]])+"\n")
