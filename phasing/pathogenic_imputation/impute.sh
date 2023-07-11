#!/bin/bash

#PBS -N impute_pathogenic
#PBS -l nodes=1:ppn=1
#PBS -l walltime=40:00:00
#PBS -o /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/phasing/pathogenic_imputation/log
#PBS -e /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/phasing/pathogenic_imputation/log
#PBS -A gymreklab-group
#PBS -q hotel
#PBS -V



source /projects/ps-gymreklab/helia/ensembl/venv_3.9/bin/activate

DIR=/projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/phasing/validation

cd /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/phasing/pathogenic_imputation/

echo $pop
echo $chr

while read sample; do
	./L1O.sh -s $sample -b "$DIR"/beagle.19Apr22.7c0.jar -c $chr -p $pop
done < "$DIR"/"$pop"_names.txt 


