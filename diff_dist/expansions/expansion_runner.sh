#!/bin/bash

#PBS -N diff_extract
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -o /projects/ps-gymreklab/helia/ensembl/experiments/charact/diff_dist/log
#PBS -e /projects/ps-gymreklab/helia/ensembl/experiments/charact/diff_dist/log
#PBS -A gymreklab-group
#PBS -q hotel
#PBS -V


source /projects/ps-gymreklab/helia/ensembl/venv_3.9/bin/activate
cd /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/diff_dist/expansions


python3 expansion_nm.py $chr
