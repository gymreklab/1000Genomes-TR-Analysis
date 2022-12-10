import pandas as pd


addr = "/projects/ps-gymreklab/helia/ensembl/experiments/coding_regions"
genes = pd.read_csv(addr + "/coding_regions.bed", delim_whitespace=True, header=None)
genes = genes[genes[0].str.len() < 6] # Main chromosomes
genes.to_csv(addr + "/coding_regions_corrected.bed", sep = "\t", header = None, index = None)
