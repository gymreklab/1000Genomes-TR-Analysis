"""Analyze kruscal-wallis results."""
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


if __name__ == '__main__':
	results_path = os.path.join('computed_stats', 'kruskal.csv')
	alpha = 1e-9

	# Load data
	results = pd.read_csv(results_path)

	# Drop STRs where all populations had the same repeat count
	results = results.dropna().reset_index(drop=True)

	# Correct p-values
	results['kruskal_pval_bonf'] = results.kruskal_pval * len(results)

	# Identify significant correlations
	results['is_significant'] = results.kruskal_pval_bonf < alpha

	# Plot ratio of signif different by chromesome
	sns.countplot(
		data=results,
		x='chr',
		hue='is_significant'
	)
	plt.show()