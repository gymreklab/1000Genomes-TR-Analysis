chr=$1
pops=('AFR' 'EUR' 'EAS' 'AMR' 'SAS')

for pop in "${pops[@]}"
do
	#qsub tag_job.sh -v n="0",pop=$pop,chr=$chr
	for i in {1..9}
	do
		qsub tag_job.sh -v n=$i,pop=$pop,chr=$chr
	done
done
