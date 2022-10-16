"""Stratify and plot differences in copy number for coding vs noncoding STRs"""

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm


if __name__ == '__main__':
	results_path = os.path.join('computed_stats', 'kruskal.csv')
	coding_data_path = os.path.join('..', 'stats', 'coding_overlap', 'TR_intersect.txt')
	alpha = 1e-9

	# Load kruskal results data
	results = pd.read_csv(results_path)

	# Drop STRs where all populations had the same repeat count
	results = results.dropna().reset_index(drop=True)

	# Correct p-values
	results['kruskal_pval_bonf'] = results.kruskal_pval * len(results)

	# Identify significant correlations
	results['Significant Difference'] = results.kruskal_pval_bonf < alpha

	# Load coding/noncoding data
	coding_data = pd.read_csv(coding_data_path, sep='\t', header=None)
	coding_data = coding_data.rename(columns={0: 'chr', 1: 'start', 2: 'end'})

	# Check if each STR is in coding region by comparing with all coding
	# regions on the same chromosome
	results['Coding'] = False
	for chr in tqdm(results.chr.unique()):
		chr_mask = results.chr == chr
		chr_coding_data = coding_data[coding_data.chr == chr]

		for i, row in tqdm(results[chr_mask].iterrows(), total=chr_mask.sum(), desc=chr):
			# Check if any coding region overlaps with this STR
			overlaps = chr_coding_data.start <= row.position
			overlaps = overlaps & (chr_coding_data.end >= row.position)
			if any(overlaps):
				results.loc[i, 'Coding'] = True

	# Save results
	results.to_csv(os.path.join('computed_stats', 'kruskal_coding.csv'), index=False)
	