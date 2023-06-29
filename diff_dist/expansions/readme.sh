#Extract bp diff from reference allele for each chromosome:

for i in {1..22}; do qsub diff_extractor.pbs -v chr=$i;done

#Analyze each locus expansions and save results:

for i in {1..22}; do qsub expansion_runner.sh -v chr=$i; done


