import pandas as pd
from collections import defaultdict
import sys

chromosome = sys.argv[1]
output_addr = "/projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/webstr/tables"

## Forming repeat tables
def fix_start(loc):
    if loc in gene_info_only:
        loc = loc.split(":")
        chrom = loc[0]
        pos = loc[1].split("-")
        start = int(pos[0]) - 1
        if len(pos) == 1:
            end = int(start) + 1
        else:
            end = pos[1]
        return chrom + ":" + str(start) + "-" + str(end)
    return loc

## Reading gene annotations
repeat_info = pd.read_csv(f"/projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/webstr/genome_annotations/gene_annotation_chr{chromosome}.txt", comment = "#", delim_whitespace=True, header=None)
repeat_info.columns = "#Uploaded_variation     Location        Allele  Gene    Feature Feature_type    Consequence     cDNA_position   CDS_position    Protein_position        Amino_acids     Codons  Existing_variation      Extra".split()
repeat_info = repeat_info.drop_duplicates(subset = ['#Uploaded_variation', 'Location', 'Gene']) # Drop transcript information

ensemble_output = pd.read_csv(f"{output_addr}/repeat_tables/ensemble_{chromosome}.txt", sep = "\t", header = None)
ensemble_output['id'] = ensemble_output[0] + ":" + ensemble_output[1].astype(str) + "-" + ensemble_output[2].astype(str)
ensemble_output = ensemble_output.drop_duplicates(subset='id')

gene_info = pd.merge(ensemble_output, repeat_info, left_on="id", right_on="Location", indicator=True, how="outer")

# Find and Fix coordinates changed by 1 by VEP
gene_info_only = list(gene_info[gene_info['_merge'] == "right_only"]['Location']) # start moved by 1
repeat_info['Location'] = repeat_info.apply(lambda x: fix_start(x['Location']), axis = 1) 


gene_info = pd.merge(ensemble_output, repeat_info, left_on="id", right_on="Location", how="outer")
gene_info = gene_info[[0,1,2,3,'id','Gene']]
gene_info.columns = ['Chrom', 'Start', 'End', 'Motif', 'ID', 'Gene']
gene_info_grouped = gene_info.groupby(['Chrom', 'Start', 'End', 
                                       'Motif', 'ID'], as_index = False).agg({'Gene':lambda x: list(x)})
gene_info_grouped['Gene'] = gene_info_grouped.apply(lambda x: json.dumps(x['Gene']), axis = 1)
gene_info_grouped['Source'] = 'EnsembleTR'


## Forming allele frequency and heterozyosity table
def fix_freqs(freqs):
    updated_freqs = defaultdict(int)
    if freqs == ".":
        return "."
    else:
        freqs = freqs.split(",")
        for freq in freqs:
            freq = freq.split(":")
            updated_freqs[len(freq[0])] += float(freq[1])
    return dict(updated_freqs)
            
pops = ['AFR', 'AMR', 'EAS', 'SAS', 'EUR']
addr = '/projects/ps-gymreklab/helia/ensembl/experiments/allele_freq/freqs/'
all_pop_df = gene_info_grouped[['Chrom', 'Start','ID']]
for pop in pops:
    freq_het = pd.read_csv(f"{addr}freqs_chr{chromosome}_{pop}.tab", sep = "\t")
    freq_het.end = freq_het.end - 1
    freq_het['ID'] = freq_het.chrom + ":" + freq_het.start.astype(str) + "-" + freq_het.end.astype(str)
    freq_het['afreq-1'] = freq_het.apply(lambda x:fix_freqs(x['afreq-1']), axis = 1)
    freq_het['afreq-1'] = freq_het.apply(lambda x: json.dumps(x['afreq-1']), axis = 1)
    freq_het = freq_het[['chrom', 'start', 'ID', 'afreq-1', f'het-1', 'numcalled-1']]
    freq_het.columns = ['Chrom', 'Start', 'ID', f'afreq_{pop}', f'het_{pop}', f'numcalled_{pop}']
    freq_het = freq_het.drop_duplicates(subset='ID')
    all_pop_df = pd.merge(all_pop_df, freq_het, on=['Chrom', 'Start','ID'], how='left')
    
    
assert(len(gene_info_grouped) == len(all_pop_df))

gene_info_grouped.to_csv(f"{output_addr}/repeat_tables/repeat_info_chr{chromosome}.csv", index=False, sep = "\t", quoting=csv.QUOTE_NONE, escapechar='', quotechar='')
all_pop_df.to_csv(f"{output_addr}/afreq_het_tables/afreq_het_chr{chromosome}.csv", index=False, sep = "\t", quoting=csv.QUOTE_NONE, escapechar='', quotechar='')
