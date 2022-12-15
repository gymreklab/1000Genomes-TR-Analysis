#!/bin/bash

# Results from Yang here: /storage/yal084/geuvidas_eQTL/results/

# Download gtex data
wget https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6917484/bin/NIHMS1540630-supplement-Sup_Dat2.gz
tar -xzvf NIHMS1540630-supplement-Sup_Dat2.gz

# Liftover the gtex coords to hg38
./liftover_gtex.sh
