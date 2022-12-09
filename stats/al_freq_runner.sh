subpops=('H3Africa' 'AFR' 'AMR' 'EUR' 'EAS' 'SAS')
for subpop in "${subpops[@]}"
do
  for chrom in {1..22}
  do
        echo $subpop $chrom
	      qsub al_freq.pbs -v pop_name="$subpop",chr=$chrom
	done
done
