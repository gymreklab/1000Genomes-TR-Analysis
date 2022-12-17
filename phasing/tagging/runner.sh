pops=('AFR' 'EUR' 'EAS' 'AMR' 'SAS')

for i in {1..9}
do
  for pop in "${pops[@]}"
  do
	qsub tag_job.sh -v n=$i,pop=$pop
  done
done
