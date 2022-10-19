"""
Plot pairwise percent different distributions between super populations
using Conover's test results.
"""

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


if __name__ == '__main__':
	stats_dir = 'computed_stats'

	alpha = 1e-9

	# Load data
	conover_data = pd.read_csv(
		os.path.join(stats_dir, 'conover.csv'), index_col=0
	)

	# Bonferroni correction
	conover_data['conover_pval_bonf'] = conover_data.conover_pval * len(conover_data)

	# Identify significant correlations
	conover_data['signif'] = conover_data.conover_pval_bonf < alpha

	# Group by pop1 and pop2 and get percent different
	percent_diff = conover_data.groupby(['pop1', 'pop2']).signif.mean().reset_index()

	# Pivot and fill diagonal with nan
	percent_diff = percent_diff.pivot(index='pop1', columns='pop2', values='signif')
	np.fill_diagonal(percent_diff.values, np.nan)

	# Plot
	sns.heatmap(percent_diff, cmap='Blues', annot=True, fmt='.2f', vmin=0)
	plt.suptitle(
		"Pairwise Differences in Copy Number at Overall Significantly Different STRs"
	)
	plt.xlabel(None)
	plt.ylabel(None)

	# Label colorbar
	cbar = plt.gcf().axes[-1]
	cbar.set_ylabel(
		r'Fraction Significantly Different ($\alpha$ = {:.1e})'.format(alpha),
		rotation=270,
		labelpad=15
	)

	plt.show()