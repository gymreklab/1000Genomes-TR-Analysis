#!/usr/bin/env python3
"""
Merge capillary and WGS calls
Learn offsets
Output table to compare

Usage:
./merge_cap_wgs.py <product_sizes> <gangstr_calls> <hipstr_calls> <ensemble_calls>
"""

import sys

try:
    psizes_file = sys.argv[1]
    gangstr_file = sys.argv[2]
    hipstr_file = sys.argv[3]
    ensemble_file = sys.argv[4]
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

from utils import *

############### Load HipSTR, GangSTR, EnsembleTR calls ###########
loci = pd.read_csv("1000g_loci.csv")
loci = loci.drop_duplicates(keep='first')
samples = [item.strip() for item in open("samples.txt", "r").readlines()]

gangstr_df = load_vcf(gangstr_file, samples, "GangSTR", gang_dict, loci)
hipstr_df = load_vcf(hipstr_file, samples, "HipSTR", hip_dict, loci)
ensemble_df = load_vcf(ensemble_file, samples, "Ensemble", ensemble_dict, loci)

############### Load capillary ###########
cap = pd.read_csv(psizes_file)
cap = pd.melt(cap, id_vars=["PrimerID","RefProductSize"], value_vars=samples, \
        var_name="sample", value_name="Cap")
cap["batch"] = cap.apply(lambda x: "PG" if x["sample"] in pg_samples else "Coriell", 1)

############### Learn offsets ###########
offsets = learn_offsets(cap, hipstr_df, gangstr_df, loci)

############### Merge everything ###########
merged = pd.merge(cap, offsets, on=["PrimerID","batch"])
hipstr_df["sample"] = hipstr_df["SampleID"]
gangstr_df["sample"] = gangstr_df["SampleID"]

merged = pd.merge(merged, hipstr_df, how="left", on=["PrimerID","SampleID"])
merged = pd.merge(merged, gangstr_df, how="left", on=["PrimerID","SampleID"])

print(merged.head())

