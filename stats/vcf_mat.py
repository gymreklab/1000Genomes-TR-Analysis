import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
import sys
n_samples = 3202


def line_to_row(line):
    row = []
    calls = line.split("\t")
    alts = calls[4]
    ref = calls[3]
    if alts == "":
        print("oops, empty alt allele at " + calls[1])
        exit()
    if alts != ".":
        alts = alts.split(",")
        alleles = [ref] + alts
    else:
        alleles = [ref]
    calls = calls[9:]
    for call in calls:
        call = call.split(":")
        gt = call[0]
        if gt == ".":
            row.append(np.nan)
        else:
            gt = gt.split("/")
            row.append(int(len(alleles[int(gt[0])]) + len(alleles[int(gt[1])])))
    assert(len(row) == n_samples)
    return row

# reading input
samples = []

with open("/projects/ps-gymreklab/helia/ensembl/experiments/pca/mat_all.txt",'w') as f:
    while True:
        line = sys.stdin.readline()
        if not line:
           break
        if line[0] == "#":
            if line[1] == "C":
                samples = line.split("\t")
                samples = samples[9:]
                assert(len(samples) == n_samples)
            else:
                continue
        else:
            row = np.array(line_to_row(line))
            col_mean = np.nanmean(row, axis=0)
            inds = np.where(np.isnan(row))
            row[inds[0]] = col_mean
            row = [str(int(x)) for x in row]
            f.write("\t".join(row))
            f.write("\n")


with open("/projects/ps-gymreklab/helia/ensembl/experiments/pca/names.txt",'w') as f2:
    f2.write(",".join(samples))





