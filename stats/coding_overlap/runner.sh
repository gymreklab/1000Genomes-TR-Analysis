for i in {1..22}
do
	qsub coding_tr_heterozygosity.sh  -v chr=$i
done
