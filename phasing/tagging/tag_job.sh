#!/bin/bash

#PBS -N tag
#PBS -l nodes=1:ppn=1
#PBS -l walltime=120:00:00
#PBS -o /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/phasing/tagging/log
#PBS -e /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/phasing/tagging/log
#PBS -A gymreklab-group
#PBS -q home
#PBS -V

cd /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/phasing/tagging
source /projects/ps-gymreklab/helia/ensembl/venv_3.9/bin/activate

echo $chr
echo $n
echo $pop

./run_tagger.sh tag_regions/chr"$chr"_0"$n" $pop $chr > files/chr"$chr"_"$n"_"$pop"_tag_info.txt
