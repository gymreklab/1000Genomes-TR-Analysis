import pandas as pd
import numpy as np
from sklearn.decomposition import IncrementalPCA
from matplotlib import pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.patches as mpatches
import sys

# mat = []
n_samples = 1241

ipca = IncrementalPCA(n_components=2)
for i in range(0,n_samples,100):
     print(i)
     mat = pd.read_csv("/projects/ps-gymreklab/helia/ensembl/experiments/pca/mat_all_AFR.txt", usecols=range(i,min(i+100,n_samples)), delimiter="\t", header = None)
     mat = mat.to_numpy()
     ipca.partial_fit(mat.T)
print("fit_finished")

x_transform = np.ndarray(shape=(0, 2))
for i in range(0,n_samples,100):
    print(i)
    mat = pd.read_csv("/projects/ps-gymreklab/helia/ensembl/experiments/pca/mat_all_AFR.txt", usecols=range(i,min(i+100,n_samples)), delimiter="\t", header = None)
    mat = mat.to_numpy()
    partial_x_transform = ipca.transform(mat.T)
    x_transform = np.vstack((x_transform, partial_x_transform))



np.savetxt("/projects/ps-gymreklab/helia/ensembl/experiments/pca/transformed_mat_AFR.txt", x_transform)
print("pca finished")


#x_transform = np.loadtxt("/projects/ps-gymreklab/helia/ensembl/experiments/pca/transformed_mat.txt")