"""Run Kruskal-Wallis H-test tests."""

import os
from collections import namedtuple

import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm


if __name__ == '__main__':
	data_dir_path = 'preprocessed_data'

	chroms = list(range(1, 23))
	# super_pops = ['AMR', 'AFR', 'EAS', 'EUR', 'SAS']

	results = []

	NullResult = namedtuple('NullResult', 'statistic pvalue')

	n_skipped = 0

	for chrom in tqdm(chroms):
		print("\nLoading CSV...")
		df = pd.read_csv(os.path.join(data_dir_path, f'chr{chrom}.csv')).drop_duplicates()
		print("\tcomplete")

		for _, row in tqdm(df.iterrows(), desc=f"chr{chrom}", total=len(df)):
			all_diffs = [
				row.diffs_AMR, row.diffs_AFR, row.diffs_EAS, 
				row.diffs_EUR, row.diffs_SAS
			]
			if any(d == '[]' for d in all_diffs):
				# print("No vals for at least one populations")
				n_skipped += 1
				continue
			all_diffs = [
				np.array(list(int(c.strip()) for c in d.strip(' []').split(','))) for d in all_diffs
			]
			
			# Run test
			try:
				kruskal_res = stats.kruskal(*all_diffs, nan_policy='raise')
			except ValueError as e:
				if str(e) == "All numbers are identical in kruskal":
					kruskal_res = NullResult(np.nan, np.nan)
				else:
					raise

			results.append({
				'chr': row.chr,
				'position': row.position,
				'motif_len': row.motif_len,
				'kruskal_statistic': kruskal_res.statistic,
				'kruskal_pval': kruskal_res.pvalue
			})

		print(f"{n_skipped} STRs skipped for missing data for at least one pop")

	# Save Results
	results = pd.DataFrame(results)
	results.to_csv(os.path.join('computed_stats', 'kruskal.csv'))

