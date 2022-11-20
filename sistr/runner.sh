for i in {1..22}
do
	qsub period_ru.pbs -v chr=$i
done
