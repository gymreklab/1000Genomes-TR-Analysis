Computing basic stats on EnsembleTR calls

VCFs at `/gymreklab-tscc/helia/ensembl/ensemble_out/`

How many calls:

wc -l `/projects/ps-gymreklab/helia/ensembl/experiments/upset_plot/methods_chr**` - 22 (#chromosomes)



### PCA

Convert the EnsembleTR vcf file for the whole genome to a matrix of #loci * #samples where each cell shows len(allele1) + len(allele2) with **vcf_mat.py**.

Performing incrementalPCA and plotting the PCA plot with **pca_plot.py**.

Submitting a job for above scripts with **pca.pbs**.


### Upset plot

Extract method per locus with **method_extract.pbs**.

Generate upset plot with **upset_plot.py**.


### Per population variants plot

Extract genotypes per sample with **variant_pop.pbs** and **qc.py**.

Making the plot with **variants_per_pop.ipynb**.

