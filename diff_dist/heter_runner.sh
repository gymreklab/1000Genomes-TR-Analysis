pops=('AFR' 'EUR')
for i in {1..22}
do
  for pop in "${pops[@]}"
  do
          qsub heter.sh -v chr="$i",pop_name="$pop"
        done
done
