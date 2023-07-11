chrom=$1
pos=$2

statSTR \
    --vcf /projects/ps-gymreklab/helia/ensembl/ensemble_out/ensemble_"$chrom"_filtered.vcf.gz \
    --vcftype hipstr \
    --afreq --use-length \
    --region "$chrom":"$pos"-"$pos" \
    --samples ../het-stats/EUR_samples.txt,../het-stats/EAS_samples.txt,../het-stats/SAS_samples.txt,../het-stats/AMR_samples.txt,../het-stats/AFR_samples.txt,../het-stats/h3africa_samples.txt \
    --sample-prefixes EUR,EAS,SAS,AMR,AFR,H3A \
    --out tmp/"$chrom"."$pos".stats
