chr=21
while read p; do
  qsub impute.pbs -v sample=$p,chr=$chr,pop=SAS
done < SAS_samples.txt

