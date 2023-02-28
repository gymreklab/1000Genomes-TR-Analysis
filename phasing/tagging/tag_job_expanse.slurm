#!/bin/bash

#SBATCH --partition=ind-shared
#SBATCH --account=ddp268
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=2G
#SBATCH --export=ALL

module load cpu gcc anaconda3
module load bcftools

cd /expanse/projects/gymreklab/helia/ensembl/1000Genomes-TR-Analysis/phasing/tagging

chr=$1
pop=$2
n=$3

./run_tagger_expanse.sh tag_regions/chr"$chr"_0"$n" $pop $chr > files/chr"$chr"_"$n"_"$pop"_tag_info.txt


#./run_tagger_expanse.sh tag_regions/chr"$chr"_0"$n"_last $pop $chr > files/chr"$chr"_"$n"_"$pop"_tag_info_last.txt
