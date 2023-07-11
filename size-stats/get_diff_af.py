#!/usr/bin/env python3

import cyvcf2
import numpy as np
import pandas as pd
import sys
import trtools.utils.tr_harmonizer

chrom = sys.argv[1]

# Load VCF
reader = cyvcf2.VCF("/projects/ps-gymreklab/helia/ensembl/ensemble_out/ensemble_chr%s_filtered.vcf.gz"%chrom)

for v in reader:
    if v.FILTER is not None: continue
    trrecord = trtools.utils.tr_harmonizer.HarmonizeRecord(trtools.utils.tr_harmonizer.VcfTypes.hipstr, v)
    afreqs = trrecord.GetAlleleFreqs(uselength=True)
    reflen = trrecord.ref_allele_length
    for a in afreqs:
        if afreqs[a] < 0.001: continue 
        if a == reflen : continue
        sys.stdout.write("\t".join([v.CHROM, str(v.POS), str(v.INFO["PERIOD"]), "%.2f"%(a-reflen), "%.3f"%afreqs[a]])+"\n")
