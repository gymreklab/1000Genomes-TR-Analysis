chrom=21
for i in regions/"$chrom"/sec_*
do
	echo $i
        qsub phase.pbs -v chr=$chrom,file=$i
done
