for i in {1..22}
do
	qsub readme.sh -v chrom=$i
done
