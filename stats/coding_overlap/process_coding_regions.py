import pandas as pd


addr = "/projects/ps-gymreklab/helia/ensembl/experiments/coding_regions"
genes = pd.read_csv(addr + "/coding_regions.bed", delim_whitespace=True)
genes = genes.drop_duplicates(subset = ['chrom', 'cdsStart', 'cdsEnd'])
genes = genes[genes['chrom'].str.len() < 6] # Main chromosomes
genes = genes[genes['cdsStart'] != genes['cdsEnd']] # Remove non-coding genes

genes = genes[["chrom", "cdsStart", "cdsEnd", "#name", "proteinID", "strand"]]
genes.to_csv(addr + "/coding_regions_corrected.bed", sep = "\t", header = None, index = None)
