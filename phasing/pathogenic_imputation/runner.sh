chroms=(4 6 14)
for pop in "AFR" "EUR" "SAS" "EAS" "AMR"
do
        for chr in "${chroms[@]}"
        do
		echo $chr
		qsub impute.sh -v chr=$chr,pop=$pop
        done
done
