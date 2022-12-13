#!/bin/bash

#PBS -N pca_plot
#PBS -l nodes=1:ppn=5
#PBS -l walltime=10:00:00
#PBS -o /projects/ps-gymreklab/helia/ensembl/experiments/pca/log
#PBS -e /projects/ps-gymreklab/helia/ensembl/experiments/pca/log
#PBS -V
#PBS -A gymreklab-group
#PBS -q hotel

source /projects/ps-gymreklab/helia/ensembl/venv_3.9/bin/activate



python3 /projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/stats/pca/pca_plot_AFR.py &&
 
echo "Done"
