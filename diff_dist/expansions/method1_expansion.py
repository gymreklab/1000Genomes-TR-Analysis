import pandas as pd
import sys
import numpy as np
chrom = sys.argv[1]

pedigree = pd.read_csv("/projects/ps-gymreklab/helia/TR_1000G/1000G.ped", delim_whitespace=True)
pedigree = pedigree[['SampleID','Superpopulation']]
samp_to_pop = pd.Series(pedigree.Superpopulation.values,index=pedigree.SampleID).to_dict()
H3Africa_names = pd.read_csv("/projects/ps-gymreklab/helia/H3Africa/names/H3A_Baylor_sample_country.txt", header=None, delim_whitespace=True)

for index,row in H3Africa_names.iterrows():
    samp_to_pop[row[0]] = "H3Africa"

header = []
with open(f"/projects/ps-gymreklab/helia/ensembl/experiments/charact/diff_dist/diff/{chrom}_diff.txt") as f:
    for line in f:
        if line.startswith("chr"):
                f.close()
                break
        header.append(line.strip())

        
header = [samp_to_pop[x] for x in header]
header = ["CHROM", "POS", "PERIOD", "MOTIF"] + header


df = pd.read_csv(f"/projects/ps-gymreklab/helia/ensembl/experiments/charact/diff_dist/diff/{chrom}_diff.txt", skiprows=3552, sep = "\t", header = None)
df = df.drop(3554, axis = 1)
df.columns = header

def flatten(l):
    return [item for sublist in l for item in sublist]

def find_pop_expansion(row):
    maximums = []
    diffs_dict = {}

    for pop in ['AMR', 'AFR', 'EAS', 'EUR', 'SAS', 'H3Africa']:
        diffs = flatten([x.split("/") for x in list(row[pop])])
        diffs = [int(int(d) / row['PERIOD']) for d in diffs if d != "."]
        if len(diffs) == 0:
            continue
        diffs_dict[pop] = diffs
        maximums.append((max(diffs), pop))

    if len(maximums) < 2:
        return
    
    maximums.sort(reverse=True)
    if maximums[0][0] < 10:
        return

    print(maximums)
    max_diffs = flatten([x.split("/") for x in list(row[maximums[0][1]])])
    max_diffs = [int(int(d) / row['PERIOD']) for d in max_diffs if d != "."]
    print(maximums[0][0], maximums[1][0])
    print(np.mean(maximums[0][0], maximums[1][0]))
    if maximums[0][0] - maximums[1][0] > 15 and \
       len([x for x in max_diffs if x > np.mean(maximums[0][0], maximums[1][0])]) > 10:
        output = [str(x) for x in [row["CHROM"], row["POS"], row["PERIOD"], row["MOTIF"], maximums[0][0], 
                  maximums[0][1], maximums[1][0]]]
        print("\t".join(output), flush=True)
        
    return
    
for index, row in df.iterrows():
    find_pop_expansion(row)
