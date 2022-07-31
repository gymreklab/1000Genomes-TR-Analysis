import pandas as pd

genes = pd.read_csv("coding_regions_sorted.bed", delim_whitespace=True, header = None)
genes.columns = ['chr','start','end', 'ID', 'x', 'strand']
genes = genes.drop_duplicates(subset = ['chr', 'start', 'end'])
genes = genes[genes['chr'].str.len() < 6]
genes.to_csv("coding_regions_sorted_corrected.bed", sep = "\t", header = None, index = None)