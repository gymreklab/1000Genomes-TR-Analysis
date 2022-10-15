subpops=('ACB' 'ASW' 'BEB' 'CDX' 'CEU' 'CHB' 'CHS' 'CLM' 'ESN' 'FIN' 'GBR' 'GIH' 'GWD' 'IBS' 'ITU' 'JPT' 'KHV' 'LWK' 'MSL' 'MXL' 'PEL' 'PJL' 'PUR' 'STU' 'TSI' 'YRI')
for subpop in "${subpops[@]}"
do
  for chrom in {1..22}
  do
        echo $subpop $chrom
	      qsub variant_pop.pbs -v pop_name="$subpop",chr=$chrom
	done
done
