#!/usr/bin/env python3
"""
Format CODIS frequencies into a table with:
chrom, start, end, freq_AFR, freq_AMR, ...

Usage:
./format_codis.py <known freqs> <freq file 1> <freq file 2>

freq files are assumed to be named $pop.tab
"""

import pandas as pd
import os
import sys
from collections import defaultdict
from functools import reduce

try:
    knownfreqs = sys.argv[1]
    freqfiles = sys.argv[2:]
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

motifs = pd.read_csv(knownfreqs, delim_whitespace=True, usecols=range(3))[["chrom","pos","motif"]]
motifs.columns = ["chrom","start","motif"]

def len_to_cn(freqs, period, is_htt=False):
    new_freqs = defaultdict(int)
    freqs = freqs.split(",")
    for freq in freqs:
        seq, num = freq.split(":")
        if is_htt:
            cn = int(seq.count("CAG")-1)
        else:
            cn = int(len(seq)/period)
        new_freqs[cn] += float(num)
    return ",".join([str(key)+":"+str(round(new_freqs[key],3)) for key in new_freqs])

dfs = []
for f in freqfiles:
    df = pd.read_csv(f, sep="\t")
    pop = os.path.basename(f).split(".")[0].split("-")[1]
    df = pd.merge(df, motifs, on=["chrom","start"])
    df["is.htt"] = False
    df.loc[df["start"]==3074877, "is.htt"] = True
    df["freq_%s"%pop] = df.apply(lambda x: len_to_cn(x["afreq-1"], len(x["motif"]), is_htt=x["is.htt"]), 1)
    dfs.append(df[["chrom", "start", "end", "freq_%s"%pop]])

ALL = reduce(lambda  left,right: pd.merge(left, right, on=['chrom','start','end'],
                                          how='outer'), dfs)
ALL.to_csv(sys.stdout, sep="\t", index=False)
