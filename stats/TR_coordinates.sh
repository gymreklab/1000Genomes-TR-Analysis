
# Get TR coordinates
rm all_TR_coordinates.bed
for i in {1..22}
do
	bcftools query  -f '%CHROM\t%POS\t%END\n' /projects/ps-gymreklab/helia/ensembl/ensemble_out/merged_chr"$i"_sorted_ver2.vcf.gz >> all_TR_coordinates.bed
done

# Intersect with genes
bedtools intersect -a all_TR_coordinates.bed -b coding_regions_sorted_corrected.bed -u > intersect.txt
