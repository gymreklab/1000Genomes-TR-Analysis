#!/usr/bin/env python3
"""
Merge capillary and WGS calls
Learn offsets
Output table to compare

Usage:
./merge_cap_wgs.py <product_sizes> <gangstr_calls> <hipstr_calls> <ensemble_calls> <asuragen calls>
"""

import sys

try:
    psizes_file = sys.argv[1]
    gangstr_file = sys.argv[2]
    hipstr_file = sys.argv[3]
    ensemble_file = sys.argv[4]
    asuragen_file = sys.argv[5]
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
        var_name="sample", value_name="Prd")
cap["batch"] = cap.apply(lambda x: "PG" if x["sample"] in pg_samples else "Coriell", 1)

############### Learn offsets ###########
offsets = learn_offsets(cap, hipstr_df, gangstr_df, loci)

############### Merge everything ###########
merged = pd.merge(cap, offsets, on=["PrimerID","batch"])
merged["SampleID"] = merged["sample"]

merged = pd.merge(merged, hipstr_df, how="left", on=["PrimerID","SampleID"])
merged = pd.merge(merged, gangstr_df, how="left", on=["PrimerID","SampleID"])
merged = pd.merge(merged, ensemble_df, how="left", on=["PrimerID","SampleID"])
merged["Cap.Binned"] = merged.apply(lambda x: GetBinnedAlleles(x["PrimerID"], x["Prd"], x["sample"]), 1)
merged["Cap"] = merged.apply(GetCap, 1)

merged["match.hipstr"] = GetMatch(merged["HipSTR"], merged["Cap.Binned"])
merged["match.gangstr"] = GetMatch(merged["GangSTR"], merged["Cap.Binned"])
merged["match.ensemble"] = GetMatch(merged["Ensemble"], merged["Cap.Binned"])

##### First output the non-asuragen ones
merged_not_asuragen = merged[~merged["PrimerID"].isin(["HTT","C9orf72"])]
printcols = ["PrimerID","SampleID","RefProductSize","period", \
        "Prd", "Cap.Binned", \
        "HipSTR","GangSTR","Ensemble",
        "match.hipstr", "match.gangstr", "match.ensemble"]
merged_not_asuragen[printcols].to_csv("TableS4-WGSvsCapillary.csv", index=False)

#### Now deal with Asurgen ones
keepcols = ["PrimerID","SampleID", "period","HipSTR","GangSTR","Ensemble"]
merged_asuragen = merged[merged["PrimerID"].isin(["HTT","C9orf72"])][keepcols].copy()

asuragen = pd.read_csv(asuragen_file)
merged_asuragen = pd.merge(merged_asuragen, asuragen, how="left", on=["PrimerID","SampleID"])

# Dummy columns
for col in ["Prd","RefProductSize"]:
    merged_asuragen[col] = np.nan

merged_asuragen["Cap.Binned"] = merged_asuragen.apply(GetBinnedAsuragen, 1)
merged_asuragen["match.hipstr"] = GetMatch(merged_asuragen["HipSTR"], merged_asuragen["Cap.Binned"])
merged_asuragen["match.gangstr"] = GetMatch(merged_asuragen["GangSTR"], merged_asuragen["Cap.Binned"])
merged_asuragen["match.ensemble"] = GetMatch(merged_asuragen["Ensemble"], merged_asuragen["Cap.Binned"])

merged_asuragen.sort_values(["PrimerID","SampleID"])[printcols].to_csv("TableS4-WGSvsCapillary-Asuragen.csv", index=False)


