
# Get TR coordinates
rm /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/all_TR_coordinates.bed
for i in {1..22}
do
	bcftools query  -f '%CHROM\t%POS\t%END\n' /projects/ps-gymreklab/helia/ensembl/ensemble_out/merged_chr"$i"_sorted_ver2.vcf.gz >> /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/all_TR_coordinates.bed
done

# Intersect with genes
bedtools intersect -a /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/all_TR_coordinates.bed -b /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/coding_regions_sorted_corrected.bed -u > /projects/ps-gymreklab/helia/ensembl/experiments/coding_regions/intersect.txt
