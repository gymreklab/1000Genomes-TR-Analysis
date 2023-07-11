# Comparison of 1000G to Capillary electrophoresis calls

```
# Update capillary calls from Google spreadsheet
# Format into product sizes matrix for Table S2
./update_capillary.py

# Extract EnsembleTR, HipSTR, and GangSTR calls
cat 1000g_loci.csv | awk -F"," '($3!="Chrom") {print $3 "\t" $4 "\t" $5}' > capillary_loci.bed
DATADIR=/projects/ps-gymreklab/helia/ensembl
./extract_wgs_calls.sh hipstr "$(ls ${DATADIR}/1000G_calls/hipstr/chrs/*.vcf.gz)"
./extract_wgs_calls.sh gangstr "$(ls ${DATADIR}/1000G_calls/gangstr/chrs/*corrected.vcf.gz)"
./extract_wgs_calls.sh ensemble "$(ls ${DATADIR}/filtered_calls/*filtered.vcf.gz)"

# Merge and learn offsets
./merge_cap_wgs.py \
   TableS2-ProductSizes.csv \
   gangstr_cap_calls.txt \
   hipstr_cap_calls.txt \
   ensemble_cap_calls.txt \
   TableS3-Asuragen.csv

# Merge the two Table S4s
cat TableS4-WGSvsCapillary.csv > TableS4-WGSvsCapillary-combined.csv
cat TableS4-WGSvsCapillary-Asuragen.csv | grep -v PrimerID >> TableS4-WGSvsCapillary-combined.csv
```
