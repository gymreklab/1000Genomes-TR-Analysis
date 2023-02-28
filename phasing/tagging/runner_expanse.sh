#!/bin/bash
pops=('AFR' 'EUR' 'EAS' 'AMR' 'SAS')

#pops=('SAS' 'AMR')
chr=$1
for pop in "${pops[@]}"
do
	for n in $(seq 0 3)
        do
		sbatch --job-name=${chr}_${n}_${pop} --output=log/tag_${pop}_${n}_${chr}_%j_last.out tag_job_expanse.slurm  "${chr}" "${pop}" "${n}"
        done
done
