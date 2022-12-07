#!/usr/bin/env python3

import pandas as pd

# Load sample info and get sample order
kgsamples = pd.read_csv("/gymreklab-tscc/helia/TR_1000G/1000G.ped", delim_whitespace=True)
h3afrsamples = pd.read_csv("/gymreklab-tscc/helia/H3Africa/names/H3A_Baylor_sample_country.txt", \
                           sep="\t", names=["SampleID","Country"])
h3afrsamples["Superpopulation"] = "H3Africa"
h3afrsamples["Population"] = "Country"
samples = pd.concat([kgsamples[["SampleID","Population","Superpopulation"]], \
                     h3afrsamples[["SampleID","Population","Superpopulation"]]])

# Sort by correct order
def GetPopOrder(spop):
    spops = ["EUR","EAS","SAS","AMR","AFR","H3Africa"]
    return "%s_%s"%(spops.index(spop), spop)

samples["Superpopulation"] = samples["Superpopulation"].apply(GetPopOrder)
samples = samples.sort_values(["Superpopulation","Population"])



