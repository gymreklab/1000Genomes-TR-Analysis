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
	results['Significant Difference'] = results.kruskal_pval_bonf < alpha

	# Plot ratio of signif different by chromesome
	sns.countplot(
		data=results,
		x='chr',
		hue='Significant Difference'
	)
	plt.show()

	results['Motif Length'] = results.motif_len.astype(str)
	results.loc[results.motif_len >= 5, 'Motif Length'] = '5+'

	sns.countplot(
		data=results,
		x='Motif Length',
		order=['1', '2', '3', '4', '5+'],
		hue='Significant Difference',
	)
	plt.suptitle("Differences in Copy Number Between Super Populations ")
	plt.show()

	# Plot STRs 

	# Find CA10
	pos_start = 49707674
	pos_end = 50237161
	ca10 = results[results.chr == 'chr17']
	ca10 = ca10[ca10.position >= pos_start]
	ca10 = ca10[ca10.position <= pos_end]

	# Plot a significant CA10 STR
	preprocessed_dir_path = 'preprocessed_data'
	chr17_df = pd.read_csv(os.path.join(preprocessed_dir_path, 'chr17.csv'))

	# Find most significant
	ca10 = ca10.sort_values('kruskal_pval_bonf').reset_index(drop=True)
	example = ca10.iloc[0]
	# example = results[results.position == 51831667].iloc[0]
	example_data = chr17_df[chr17_df.position == example.position]

	# Get diffs for each population and make dataframe
	diff_data = []
	for pop in ['AMR', 'AFR', 'EAS', 'EUR', 'SAS']:
		diffs = example_data[f'diffs_{pop}'].item().strip(' []').split(',')
		diff_data.extend(
			[{'pop': pop, 'diff_n_bases': int(d.strip())} for d in diffs]
		)

	diff_data = pd.DataFrame(diff_data)
	diff_data['Copy Num. From Ref.'] = diff_data.diff_n_bases / example_data.motif_len.item()

	# Plot
	sns.displot(
		data=diff_data,
		x='Copy Num. From Ref.',
		hue='pop',
		kind="kde",
		common_norm=False,
		common_grid=True,
		bw_adjust=1.5#diff_data['Copy Num. From Ref.'].std()
	)
	# sns.displot(
	# 	data=diff_data,
	# 	x='diff',
	# 	hue='pop',
	# 	kind="hist",
	# 	common_norm=True,
	# 	stat='percent',
	# 	discrete=True,
	# 	multiple='stack'
	# )
	plt.title("KW stat: {:.2f}    corrected p-value: {:.2e}".format(
		example.kruskal_statistic, example.kruskal_pval_bonf
	))
	plt.suptitle(f"Copy Number Differences for CA10 STR at chr17 {example.position}")
	plt.tight_layout()
	plt.show()

