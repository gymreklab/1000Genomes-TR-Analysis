chr=$1
pop=$2
while read p; do
  qsub impute.pbs -v sample=$p,chr=$chr,pop=$pop
done < "$pop"_names.txt

