#!/bin/bash


#PBS -N method1_expansion
#PBS -l nodes=1:ppn=1
#PBS -l walltime=20:00:00
#PBS -o /projects/ps-gymreklab/helia/ensembl/experiments/charact/diff_dist/log
#PBS -e /projects/ps-gymreklab/helia/ensembl/experiments/charact/diff_dist/log
#PBS -A gymreklab-group
#PBS -q hotel
#PBS -V



source /projects/ps-gymreklab/helia/ensembl/venv_3.9/bin/activate

python3 /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/diff_dist/expansions/method1_expansion.py $chr > /projects/ps-gymreklab/helia/ensembl/experiments/charact/diff_dist/diff/expansions_"$chr".txt
