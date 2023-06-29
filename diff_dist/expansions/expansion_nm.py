import pandas as pd
import numpy as np
from tqdm import tqdm
tqdm.pandas(desc='My bar!')
from collections import defaultdict
import sys


chrom = sys.argv[1]



##### Reading sample names and their populations #####
pedigree = pd.read_csv("/projects/ps-gymreklab/helia/TR_1000G/1000G.ped", delim_whitespace=True)
pedigree = pedigree[['SampleID','Superpopulation']]
samp_to_pop = pd.Series(pedigree.Superpopulation.values,index=pedigree.SampleID).to_dict()
H3Africa_names = pd.read_csv("/projects/ps-gymreklab/helia/H3Africa/names/H3A_Baylor_sample_country.txt", header=None, delim_whitespace=True)

for index,row in H3Africa_names.iterrows():
    samp_to_pop[row[0]] = "H3Africa"

for samp in samp_to_pop:
    if samp_to_pop[samp] == "H3Africa" or samp_to_pop[samp] == "AFR":
        samp_to_pop[samp] = "AFR"
    else:
        samp_to_pop[samp] = "NON_AFR"




header = []
with open(f"/projects/ps-gymreklab/helia/ensembl/experiments/charact/diff_dist/diff/{chrom}_diff.txt") as f:
    for line in f:
        if line.startswith("chr"):
                f.close()
                break
        header.append(line.strip())


header_pop = [samp_to_pop[x] for x in header]
header_pop = ["CHROM", "POS", "PERIOD", "MOTIF"] + header_pop


df = pd.read_csv(f"/projects/ps-gymreklab/helia/ensembl/experiments/charact/diff_dist/diff/{chrom}_diff.txt", skiprows=3552, sep = "\t", header = None)
df = df.drop(3554, axis = 1)
df.columns = header_pop


pops = ['AFR', 'NON_AFR']

def find_outlier(row):
    cn_dif = defaultdict(list)
    pop_n_outliar = []
    for pop in pops:
        pop_data = list(row[pop])
        for gbs in pop_data:
            if gbs != ".":
                cn_dif[pop].extend([int(x)/int(row['PERIOD']) for x in gbs.split("/")])

    ### Define outlier
    all_cns = cn_dif[pops[0]] + cn_dif[pops[1]]
    Q1 = np.percentile(all_cns, 25, interpolation = 'midpoint')
    Q3 = np.percentile(all_cns, 75,  interpolation = 'midpoint')
    IQR = Q3 - Q1
    strong_outlier_threshold = Q3 +3*IQR

    for pop in pops:
        pop_n_outliar.append(len([x for x in cn_dif[pop] if x > strong_outlier_threshold]) / len(cn_dif[pop]))

    return strong_outlier_threshold, pop_n_outliar[0], pop_n_outliar[1]

out = df.progress_apply(lambda row: find_outlier(row), axis = 1)


selected = df[['CHROM','POS','PERIOD','MOTIF']]
selected[['outlier_threshold', 'AFR_freq', 'NON_AFR_freq']] = pd.DataFrame(list(out), columns = ['x','y','z'])
selected = selected[((selected['AFR_freq'] > 0.01) | (selected['NON_AFR_freq'] > 0.01)) & (selected['outlier_threshold'] > 4)]
#selected = selected[(selected['AFR_freq'] > 10 * selected['NON_AFR_freq']) | (selected['NON_AFR_freq'] > 10 * selected['AFR_freq'])]

selected.to_csv(f"{chrom}_expansions.csv", index=False)

