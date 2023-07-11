import sys
import time
import datetime
from collections import defaultdict



pop = sys.argv[1]
chrom = sys.argv[2]

samples = []
ref_samples = defaultdict(int)
non_ref_samples = defaultdict(int)



alt_diff = defaultdict(str)
while True:
    line = sys.stdin.readline()
    if not line:
       break
    if line[0] == "#":
        if line[1] == "C":
            samples = line.split("\t")
            samples = samples[9:]
            samples = [sample.strip() for sample in samples]
        else:
            continue
    else:
        calls = line.split("\t")
        alts = calls[4]
        ref = calls[3]
        pos = calls[1]
        period = int(calls[7].split(";")[2].replace("PERIOD=", ""))
        if period == 1:
             continue
        calls = calls[9:]
        assert (len(calls) == len(samples))
        for i in range(len(calls)):
            call = calls[i].split(":")
            gt = call[0]
            if gt == ".":
                continue
            gt = gt.split("/")
            for al in gt:
                if al == "0":
                    ref_samples[samples[i]] += 1
                else:
                    non_ref_samples[samples[i]] += 1

print(ref_samples)
with open(f"/projects/ps-gymreklab/helia/ensembl/experiments/stats/files/{pop}_{chrom}_no_homo.txt",'w') as f:
    for key in ref_samples:
        f.write(f"{key}\t{ref_samples[key]}\t{non_ref_samples[key]}\n")


