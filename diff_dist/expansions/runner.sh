pops=('AFR' 'EUR' 'EAS' 'AMR' 'SAS' 'H3Africa')

for i in {1..22}
do
  for pop in "${pops[@]}"
  do
          qsub diff_extractor.pbs -v pop=$pop,chr=$i
        done
done

