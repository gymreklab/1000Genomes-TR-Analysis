"""
Compute Conover–Iman statisitics and p-vals for STRs with significant
Kruskal-Wallis test results.
"""

import os

import numpy as np
import pandas as pd
from scikit_posthocs import posthoc_conover
from tqdm import tqdm

import warnings
warnings.filterwarnings("error")


if __name__ == '__main__':
	corr_res_dir = 'computed_stats'
	data_dir_path = 'preprocessed_data'

	alpha = 1e-9
	
	# Load Kruskal-Wallis results
	kruskal_res = pd.read_csv(
		os.path.join(corr_res_dir, 'kruskal.csv'),
		index_col=0
	)
	kruskal_res = kruskal_res.dropna().reset_index(drop=True)
	kruskal_res['kruskal_pval_bonf'] = kruskal_res.kruskal_pval * len(kruskal_res)
	signif_STRs = kruskal_res[kruskal_res.kruskal_pval_bonf < alpha]

	# For testing, subset to chr 17
	# signif_STRs = signif_STRs[signif_STRs.chr == 'chr17']
	del kruskal_res

	# Load preprocessed data for each chromosome and run tests for
	# significant STRs
	pairwise_res = []

	for chrom in tqdm(signif_STRs.chr.unique(), desc='Chromosomes'):
		# Load data
		chr_df = pd.read_csv(os.path.join(data_dir_path, f'{chrom}.csv'))

		# Remoce rows with '[]' for a pop diff
		chr_df = chr_df[~chr_df.diffs_AMR.str.contains('\[\]')].reset_index(drop=True)
		chr_df = chr_df[~chr_df.diffs_AFR.str.contains('\[\]')].reset_index(drop=True)
		chr_df = chr_df[~chr_df.diffs_EAS.str.contains('\[\]')].reset_index(drop=True)
		chr_df = chr_df[~chr_df.diffs_EUR.str.contains('\[\]')].reset_index(drop=True)
		chr_df = chr_df[~chr_df.diffs_SAS.str.contains('\[\]')].reset_index(drop=True)

		# Remove duplicate rows in terms of chr and postition
		chr_df = chr_df[chr_df.position.isin(signif_STRs.position)]
		chr_df = chr_df.drop_duplicates(subset=['chr', 'position']).reset_index(drop=True)

		signif_chr_df = signif_STRs.merge(chr_df, on=['chr', 'position'], how='inner')

		# Compute pairwise Conover-Iman statistics
		for _, row in tqdm(
			signif_chr_df.iterrows(), desc="Computing pairwise differences", total=len(chr_df)
		):
			all_diffs = {
				'AMR': row.diffs_AMR, 'AFR': row.diffs_AFR, 'EAS': row.diffs_EAS,
				'EUR': row.diffs_EUR, 'SAS': row.diffs_SAS
			}
			if any(d == '[]' for d in all_diffs.values()):
				# print("No vals for at least one populations")
				continue
			all_diffs = {
				pop: list(int(c.strip()) for c in d.strip(' []').split(','))
				for pop, d in all_diffs.items()
			}
			all_diffs = pd.DataFrame([
				{'pop': pop, 'copy_num_diff': diff} for pop, diffs in all_diffs.items() for diff in diffs
			])

			# Compute pairwise Conover-Iman statistics
			pairwise_diffs = posthoc_conover(
				all_diffs,
				val_col='copy_num_diff',
				group_col='pop'
			)
			# Convert to long format
			pairwise_diffs = pairwise_diffs.stack().reset_index()
			pairwise_diffs.columns = ['pop1', 'pop2', 'conover_pval']
			pairwise_res.append(pairwise_diffs)

	# Save results
	pairwise_res = pd.concat(pairwise_res)
	pairwise_res.to_csv(os.path.join(corr_res_dir, 'conover.csv'))
			