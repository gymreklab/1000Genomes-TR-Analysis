import pandas as pd
import numpy as np
from sklearn.decomposition import IncrementalPCA
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import sys

mat = []


ipca = IncrementalPCA(n_components=2)
for i in range(0,3202,100):
    print(i)
    mat = pd.read_csv("/projects/ps-gymreklab/helia/ensembl/experiments/pca/mat_all.txt", usecols=range(i,min(i+100,3202)), delimiter="\t", header = None)
    mat = mat.to_numpy()
    ipca.partial_fit(mat.T)

print("fit_finished")

x_transform = np.ndarray(shape=(0, 2))
for i in range(0,3202,100):
    print(i)
    mat = pd.read_csv("/projects/ps-gymreklab/helia/ensembl/experiments/pca/mat_all.txt", usecols=range(i,min(i+100,3202)), delimiter="\t", header = None)
    mat = mat.to_numpy()
    partial_x_transform = ipca.transform(mat.T)
    x_transform = np.vstack((x_transform, partial_x_transform))



np.savetxt("/projects/ps-gymreklab/helia/ensembl/experiments/pca/transformed_mat.txt", x_transform)
print("pca finished")


#x_transform = np.loadtxt("/projects/ps-gymreklab/helia/ensembl/experiments/pca/transformed_mat.txt")

samples = []
with open("/projects/ps-gymreklab/helia/ensembl/experiments/pca/names.txt") as f:
    for line in f:
        samples = line.strip()
        samples = line.split(",")

pop_to_color = {"None": "gray"}
for pop in ["ACB","ASW"]: pop_to_color[pop] = "darkorange" # Admixed African
for pop in ["GIH","BEB","ITU"]: pop_to_color[pop] = "mediumpurple" # South Asian
for pop in ["CDX","CHB","CHS","JPT","KHV"]: pop_to_color[pop] = "yellowgreen" # East Asian
for pop in ["CEU","FIN","GBR","IBS","TSI"]: pop_to_color[pop] = "mediumturquoise" # European
for pop in ["CLM","MXL","PEL","PJL","PUR","STU"]: pop_to_color[pop] = "brown" # American
for pop in ["LWK","MSL","YRI","ESN","GWD"]: pop_to_color[pop] = "gold" # African

pedigree = pd.read_csv("/projects/ps-gymreklab/helia/TR_1000G/1000G.ped", delim_whitespace=True)
pedigree = pedigree[['SampleID','Population']]
samp_to_pop = pd.Series(pedigree.Population.values,index=pedigree.SampleID).to_dict()

colors = []
for sample in samples:
    colors.append(pop_to_color[samp_to_pop[sample.strip()]])

plt.figure(figsize=(10,10))
plt.scatter(x_transform[:,0], x_transform[:,1], c=colors)

blue_patch = mpatches.Patch(color='darkorange', label='Admixed African')
orange_patch = mpatches.Patch(color='mediumpurple', label='South Asian')
green_patch = mpatches.Patch(color='yellowgreen', label='East Asian')
yellow_patch = mpatches.Patch(color='mediumturquoise', label='European')
purple_patch = mpatches.Patch(color='brown', label='American')
red_patch = mpatches.Patch(color='gold', label='African')

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(handles=[blue_patch,orange_patch,green_patch, yellow_patch,purple_patch,red_patch], loc="lower left", fontsize=14)
plt.xlabel('PC 1 (%.2f%%)' % (ipca.explained_variance_ratio_[0]*100))
plt.ylabel('PC 2 (%.2f%%)' % (ipca.explained_variance_ratio_[1]*100))
# plt.xlabel('PC 1', fontsize=18, labelpad=20)
# plt.ylabel('PC 2',fontsize=18)

plt.savefig("/projects/ps-gymreklab/helia/ensembl/experiments/pca/pca.png",dpi=1200)
