### TR coordinates and intersect with coding regions

Coding regions were downloaded from UCSC tables and were filtered with **process_coding_regions.py** to remove the followings:

1. Drop duplicated regions
2. Genes located in chromosomes other than chr1-chr22

Get coordinates of all TR regions in EnsembleTR output with **TR_coordinates.sh**.

Intersect with coding regions using bedtools with **TR_coordinates.sh** 

Finding TRs entirely in coding regions with overlap_finder.cpp. Compare the results with overlaps found by bedtools in **coding_region_intersect.ipynb**

### Heterozygosity of coding TRs in different populations

Extract coding calls from vcf files with **extract_coding.sh**.

Calculate heterozygosity per population with **coding_tr_heterozygosity.sh** and **coding_tr_heterozygosity_runner.sh**

Generating supplementary table table 7 in **coding_region_intersect.ipynb**

### Polymorphism percentage

Get CHROM, POS, PERIOD, ALT of coding regions and whole genome with **wg_alt.sh** script. 

Calculate polymorphism rate in coding and whole genome in **coding_region_intersect.ipynb**


### Motif length distribution 


Get CHROM, POS, PERIOD, ALT of coding regions and whole genome with **wg_alt.sh** script. 

Generating Supplementary Fig. 6a in **coding_region_period_dist.ipynb**.


### Per sample stats


Average of TRs in each sample for which one or both alleles did not match the reference genome:

Filter non-ref calls for each sample using **whole_genome_nr.sh**

Count the numbers for whole genome and coding regions in **coding_non_ref.ipynb**.

Generating Supplementary Table 6 in **coding_non_ref.ipynb**.


