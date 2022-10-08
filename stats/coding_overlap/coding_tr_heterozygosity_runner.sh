
pops=('AFR' 'EUR' 'EAS' 'AMR' 'SAS')
for i in {1..22}
do
  for pop in "${pops[@]}"
  do
	  qsub coding_tr_heterozygosity.sh -v chr="$i",pop_name="$pop"
	done
done
