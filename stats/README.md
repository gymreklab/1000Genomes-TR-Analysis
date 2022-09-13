Computing basic stats on EnsembleTR calls

VCFs at `/gymreklab-tscc/helia/ensembl/ensemble_out/`

```
# How many unique loci were called?
for vcf in $(ls /gymreklab-tscc/helia/ensembl/ensemble_out/*.vcf.gz)
do 
    echo $vcf $(zcat $vcf | grep -v "^#" | wc -l)
done
/gymreklab-tscc/helia/ensembl/ensemble_out/merged_chr10_sorted_ver2.vcf.gz 86098
/gymreklab-tscc/helia/ensembl/ensemble_out/merged_chr11_sorted_ver2.vcf.gz 83155
/gymreklab-tscc/helia/ensembl/ensemble_out/merged_chr12_sorted_ver2.vcf.gz 93247
/gymreklab-tscc/helia/ensembl/ensemble_out/merged_chr13_sorted_ver2.vcf.gz 58603
/gymreklab-tscc/helia/ensembl/ensemble_out/merged_chr14_sorted_ver2.vcf.gz 59520
/gymreklab-tscc/helia/ensembl/ensemble_out/merged_chr15_sorted_ver2.vcf.gz 53028
/gymreklab-tscc/helia/ensembl/ensemble_out/merged_chr16_sorted_ver2.vcf.gz 58705
/gymreklab-tscc/helia/ensembl/ensemble_out/merged_chr17_sorted_ver2.vcf.gz 63200
/gymreklab-tscc/helia/ensembl/ensemble_out/merged_chr18_sorted_ver2.vcf.gz 47091
/gymreklab-tscc/helia/ensembl/ensemble_out/merged_chr19_sorted_ver2.vcf.gz 55512
/gymreklab-tscc/helia/ensembl/ensemble_out/merged_chr1_sorted_ver2.vcf.gz 153274
/gymreklab-tscc/helia/ensembl/ensemble_out/merged_chr20_sorted_ver2.vcf.gz 43107
/gymreklab-tscc/helia/ensembl/ensemble_out/merged_chr21_sorted_ver2.vcf.gz 21441
/gymreklab-tscc/helia/ensembl/ensemble_out/merged_chr22_sorted_ver2.vcf.gz 26452
/gymreklab-tscc/helia/ensembl/ensemble_out/merged_chr2_sorted_ver2.vcf.gz 151133
/gymreklab-tscc/helia/ensembl/ensemble_out/merged_chr3_sorted_ver2.vcf.gz 125191
/gymreklab-tscc/helia/ensembl/ensemble_out/merged_chr4_sorted_ver2.vcf.gz 115218
/gymreklab-tscc/helia/ensembl/ensemble_out/merged_chr5_sorted_ver2.vcf.gz 110036
/gymreklab-tscc/helia/ensembl/ensemble_out/merged_chr6_sorted_ver2.vcf.gz 108567
/gymreklab-tscc/helia/ensembl/ensemble_out/merged_chr7_sorted_ver2.vcf.gz 98068
/gymreklab-tscc/helia/ensembl/ensemble_out/merged_chr8_sorted_ver2.vcf.gz 91595
/gymreklab-tscc/helia/ensembl/ensemble_out/merged_chr9_sorted_ver2.vcf.gz 73749

sums to 1,775,990

# How many called by a single method?
```



### Polymorphism percentage and period distribution 

Get CHROM, POS, PERIOD, ALT of coding regions and whole genome with **stats.sh** script. 

Calculating (#loci with at least one alt allele / #all loci) in **coding_region_period_dist.ipynb**.

Generating Supplementary Fig 6a in **coding_region_period_dist.ipynb**.
