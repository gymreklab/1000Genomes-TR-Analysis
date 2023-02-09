#!/bin/bash

#PBS -N tag
#PBS -l nodes=1:ppn=1
#PBS -l walltime=20:00:00
#PBS -o /projects/ps-gymreklab/helia/ensembl/phasing_st/phased/log
#PBS -e /projects/ps-gymreklab/helia/ensembl/phasing_st/phased/log
#PBS -A gymreklab-group
#PBS -q hotel
#PBS -V

cd /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/phasing/tagging
source /projects/ps-gymreklab/helia/ensembl/venv_3.9/bin/activate


./run_tagger.sh tag_regions/chr"$chr"_0"$n".txt $pop $chr > files/chr"$chr"_"$n"_"$pop"_tag_info.txt
