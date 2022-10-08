for i in {1..22}
do
	qsub whole_genome_nr.sh   -v chr=$i
done
