for i in {1..22}
do
	qsub stats.sh -v chr=$i
done
