import sys
import numpy as np

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

### Reading STR information ###
str_addr = sys.argv[1]
snp_addr = sys.argv[2]
str_vector = ""
with open(str_addr) as f:
    for vector in f:
        if vector.strip().split("\t")[2] == ".":
            str_vector = vector.strip()
            f.close()
            break
if str_vector == "":
    sys.exit(0)
snps = []
with open(snp_addr) as f:
    snps = f.readlines()


str_vector = str_vector.split("\t")
str_alleles = [str_vector[3]] + str_vector[4].split(",")
str_vector_reform = []

for gt in str_vector[9:]:
    gt = gt.split(":")[0]
    if gt == ".":
        str_vector_reform.extend([np.nan, np.nan])
    else:
        gt = gt.split("|")
        gt = [int(x) for x in gt]
        str_vector_reform.extend([len(str_alleles[gt[0]]), len(str_alleles[gt[1]])])

if len(set(str_vector_reform)) == 1:
    eprint("Non polymorphic locus at " + str_vector[0] + "-" + str_vector[1])
    sys.exit(0)

LDs = []
cnt = 0
for i in range(len(snps)):
    snp_vector = snps[i]
    snp_vector = snp_vector.split("\t")
    snp_vector_reform = []
    for gt in snp_vector[9:]:
        gt = gt.split(":")[0]
        if gt == ".":
            snp_vector_reform.extend([np.nan, np.nan])
        if "|" in gt:
            gt = gt.split("|")
            gt = [int(x) for x in gt]
            snp_vector_reform.extend([gt[0], gt[1]])

    assert (len(str_vector_reform) == len(snp_vector_reform) == 3202*2)
    LDs.append((np.corrcoef(str_vector_reform, snp_vector_reform)[1, 0], snp_vector[1]))


assert(len(LDs) == len(snps))
max_snp = max(LDs,key=lambda item:abs(item[0]))
output = [str(x) for x in [str_vector[0], str_vector[1], max_snp[1], max_snp[0]]]

print("\t".join(output), flush = True)
