for i in {1..22}
do
	qsub method_extract.pbs -v chr=$i
done
