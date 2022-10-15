import os
from itertools import product

import pandas as pd
from tqdm import tqdm

if __name__ == '__main__':
	data_dir_path = '/projects/ps-gymreklab/helia/ensembl/experiments/charact/diff_dist/diff'

	chroms = list(range(1, 23))
	super_pops = ['AMR', 'AFR', 'EAS', 'EUR', 'SAS']
	
	# for chrom, pop in tqdm(
	# 	product(chroms, super_pops), 
	# 	total=len(chroms) * len(super_pops),
	# 	desc="Loading all chromosomes for all populations"
	# ):
	for chrom in tqdm(chroms, desc="For each chromosome"):
		chrom_data = {p: [] for p in super_pops}

		for pop in tqdm(super_pops, desc=f"chr{chrom} for each population"):
			with open(os.path.join(data_dir_path, f'{pop}_{chrom}_diff.txt'), 'r') as f:
				file_data = []

				while (line := f.readline().rstrip()):
					# First three positions of line are chr, pos, and motif len
					line = line.split('\t')
					chrom_str = line.pop(0)
					position = line.pop(0)
					motif_len = line.pop(0)

					# Pool copy number difference counts
					diff_counts = []
					for sample in line:
						sample = sample.split(':')[1]
						if sample == '.':
							continue
						else:
							diff_counts.extend(sample.split('/'))

					file_data.append({
						'chr': chrom_str,
						'position': position,
						'motif_len': motif_len,
						f'diffs_{pop}': [int(c) for c in diff_counts]
					})
		
			chrom_data[pop].append(pd.DataFrame(file_data))

		# Reformat
		chrom_data = {k: pd.concat(chrom_data[k], ignore_index=True) for k in chrom_data.keys()}

		merged_data = chrom_data.pop(super_pops[0])
		for pop in chrom_data.keys():
			merged_data = merged_data.merge(
				chrom_data[pop],
				on=['chr', 'position', 'motif_len'],
				how='outer'
			)

		# Save
		merged_data.to_csv(
			os.path.join('preprocessed_data', f'chr{chrom}.csv'),
			index=False
		)

