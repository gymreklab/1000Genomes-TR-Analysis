import pandas as pd
import sys

freq_addr = "/projects/ps-gymreklab/helia/ensembl/experiments/allele_freq/freqs"
addr="/projects/ps-gymreklab/helia/ensembl/sistr_analysis"

pop = sys.argv[1]
all_chrom = pd.DataFrame(columns = ['chrom', 'start', 'end', 'acount-1', 'numcalled-1', 'period', 'motif'])
for chrom in range(1,23):
    ru_motif = pd.read_csv(f"{addr}/al_freqs/info_{chrom}.txt", header=None, sep="\t")
    ru_motif.columns = ['chrom', 'start', 'period', 'motif']
    ru_motif = ru_motif.drop_duplicates(subset=['chrom', 'start'], keep="first")


    freqs = pd.read_csv(f"{freq_addr}/freqs_chr{chrom}_{pop}.tab", sep="\t")
    freqs = freqs.drop_duplicates(subset=['chrom', 'start'], keep="first")


    all = pd.merge(ru_motif, freqs, on=['chrom', 'start'], how='inner')
    all['numcalled-1'] = pd.to_numeric(all['numcalled-1'])
    all['numcalled-1'] = all['numcalled-1'].apply(lambda x: x*2)
    all = all[['chrom', 'start', 'end', 'acount-1', 'numcalled-1', 'period', 'motif']]
    all = all[all['acount-1'] != "."]
    all = all[(all['period'] != 1) & (all['period'] < 5)]
    all_chrom = pd.concat([all_chrom, all])

all_chrom.to_csv(f"{addr}/al_freqs/{pop}_sistr_input.txt", header=False, index=False, sep="\t")
