for i in {1..22}
do
	qsub script.sh -v chr=$i
done
