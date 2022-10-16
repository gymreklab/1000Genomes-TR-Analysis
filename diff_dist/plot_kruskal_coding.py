"""Plot signficant differences between populations for coding vs noncoding STRs"""

import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
import seaborn as sns


if __name__ == '__main__':
	results_path = os.path.join('computed_stats', 'kruskal_coding.csv')
	alpha = 1e-9

	# Load data
	results = pd.read_csv(results_path)
	results['Motif Length'] = results.motif_len.astype(str)
	results.loc[results.motif_len >= 5, 'Motif Length'] = '5+'

	# Identify significant correlations
	results['Significant Difference'] = results.kruskal_pval_bonf < alpha

	# Plot
	sns.catplot(
		kind='count',
		data=results,
		x='Motif Length',
		order=['1', '2', '3', '4', '5+'],
		hue='Significant Difference',
		col='Coding',
		sharey=False,
	)
	plt.suptitle("Differences in Copy Number Between Super Populations at Heterozygous STRs")
	plt.tight_layout()
	plt.show()

