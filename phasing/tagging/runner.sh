chr=$1
pops=('AFR' 'EUR' 'EAS' 'AMR' 'SAS')

for pop in "${pops[@]}"
do
	for i in $(seq 0 2)
	do
		qsub tag_job.sh -v n=$i,pop=$pop,chr=$chr
	done
done
