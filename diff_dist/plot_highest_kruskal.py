"""
Plot STR copy number difference distributions for ones with highest
Kruskal statistic.
"""

import os

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import trange


if __name__ == '__main__':
	results_path = os.path.join('computed_stats', 'kruskal.csv')
	data_dir_path = 'preprocessed_data'
	plot_save_dir = 'top_kruskal_plots'

	plot_top_n = 250

	# Load data
	results = pd.read_csv(results_path, index_col=0)

	# Drop STRs where all populations had the same repeat count
	results = results.dropna().reset_index(drop=True)

	# Correct p-values
	results['kruskal_pval_bonf'] = results.kruskal_pval * len(results)

	# Sort by kruskal statistic
	results = results.sort_values('kruskal_statistic', ascending=False).reset_index(drop=True)

	# For testing only chr17
	results = results[results.chr == 'chr17']

	str_data_cache = {}

	# Plot top n
	for i in trange(plot_top_n):
		# Get STR data
		if results.iloc[i].chr not in str_data_cache:
			str_data = pd.read_csv(os.path.join(data_dir_path, f'{results.iloc[i].chr}.csv'))
			str_data_cache[results.iloc[i].chr] = str_data.drop_duplicates()
		str_data = str_data_cache[results.iloc[i].chr]

		example_data = str_data[str_data.position == results.iloc[i].position].iloc[0]

		# Get diffs for each population and make dataframe
		diff_data = []
		for pop in ['AMR', 'AFR', 'EAS', 'EUR', 'SAS']:
			diffs = example_data[f'diffs_{pop}'].strip(' []').split(',')
			diff_data.extend(
				[{'Super Population': pop, 'diff_n_bases': int(d.strip())} for d in diffs]
			)

		diff_data = pd.DataFrame(diff_data)
		diff_data['Copy Num. From Ref.'] = diff_data.diff_n_bases / example_data.motif_len

		# Plot
		sns.displot(
			data=diff_data,
			x='Copy Num. From Ref.',
			hue='Super Population',
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
			results.iloc[i].kruskal_statistic, results.iloc[i].kruskal_pval_bonf
		))
		plt.suptitle(f"STR Copy Number Distribution for {results.iloc[i].chr}:{results.iloc[i].position}")
		plt.tight_layout()
		plt.savefig(os.path.join(plot_save_dir, f"{results.iloc[i].chr}_{results.iloc[i].position}.png"))
		plt.close()
