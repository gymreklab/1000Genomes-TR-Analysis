import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

pop = argv[1]
chr = argv[2]
samples = []

pop_dict = {"EUR":"European", "AFR": "African"}
pop_color = {'EUR':"blue", "AFR":"orange"}

with open(pop + "_names.txt") as f:
    for line in f:
        samples.append(line.strip())


vcf_addr=f"/projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/phasing/snps/phased_strs_chr{chr}_normalized.vcf.gz"

# Load all POS from original vcf file
all_imputed = pd.read_csv(vcf_addr, header=None,delim_whitespace=True, usecols=[1],comment='#')
all_gt = pd.read_csv(vcf_addr, header=None,delim_whitespace=True, usecols=[1],comment='#')
all_imputed.columns = ["pos"]
all_gt.columns = ["pos"]
all_imputed['pos'] = all_imputed['pos'].astype('int').astype('str')
all_gt['pos'] = all_gt['pos'].astype('int').astype('str')
all_gt = all_gt.drop_duplicates(subset='pos', keep="first")
all_imputed = all_imputed.drop_duplicates(subset='pos', keep="first")

for sample in samples:
    print(sample)
    df = pd.read_csv("diff_" + pop + "/" + chr + "/" + sample + "." + chr + ".diff.txt", header = None, delim_whitespace=True)
    df.columns = ['pos', 'gt1_' + sample, 'gt2_' + sample, 'i1_' + sample, 'i2_' + sample]
    df['gt_' + sample] = list(zip(df['gt1_' + sample], df['gt2_' + sample]))
    df['i_' + sample] = list(zip(df['i1_' + sample], df['i2_' + sample]))
    df['pos'] = df['pos'].astype('int').astype('str')
    df = df.drop_duplicates(subset='pos', keep="first")
    all_imputed = pd.merge(all_imputed, df[['pos', 'i_' + sample]], on='pos', how='left')
    all_gt = pd.merge(all_gt, df[['pos', 'gt_' + sample]], on='pos', how='left')

all_gt = all_gt.dropna(axis = 0, thresh=20)
all_imputed = all_imputed.dropna(axis = 0, thresh=20)

all_gt.to_csv("diff_" + pop + "/" + chr + "/" + "genotypes.csv", index=False)
all_imputed.to_csv("diff_" + pop + "/" + chr + "/" + "imputed.csv", index=False)


