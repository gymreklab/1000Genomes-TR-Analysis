#!/usr/bin/env python3

import numpy as np
import pandas as pd
import gspread
from oauth2client.service_account import ServiceAccountCredentials
scope = ['https://spreadsheets.google.com/feeds',
         'https://www.googleapis.com/auth/drive']

credentials = ServiceAccountCredentials.from_json_keyfile_name(
    '/storage/mgymrek/workspace/capillary-electrophoresis-webserver/capillaryelectrophoresis-e683d5dd37a2.json', scope)
gc = gspread.authorize(credentials)

def LoadGSheet(wks):
    data = wks.get_all_values()
    headers = data.pop(0)
    return pd.DataFrame(data, columns=headers)

############### Load data from the GSheet ###########
sampdata = LoadGSheet(gc.open("1000GenomesRepeatValidationDatabase").worksheet("samples"))

############### Load list of loci to consider ###########
loci = pd.read_csv('1000g_loci.csv')
primer_ids = set(loci["LocusID"])

############### Load product sizes from google sheet ###########
psdata = LoadGSheet(gc.open("1000GenomesRepeatValidationDatabase").worksheet("product size"))
psdata = psdata[psdata["PrimerID"].isin(primer_ids)][["PrimerID","SampleID","Product size"]]
psdata.columns = ["PrimerID", "SampleID","ProductSize"]
psdata["SampleID"] = psdata["SampleID"].apply(lambda x: x.strip())

# Sample product sizes
samp_prod_sizes = psdata[~psdata["SampleID"].isin(["Reference", "reference"])]
samp_prod_sizes = samp_prod_sizes[samp_prod_sizes["ProductSize"].apply(lambda x: "/" in x and x.split("/")[1].strip() != "")]
samp_prod_sizes[["prd_size_1", "prd_size_2"]] = samp_prod_sizes["ProductSize"].str.split('/', expand=True)[[0,1]].astype(float).astype(int)

############### Load ref product sizes ###########
ref_prod_sizes = LoadGSheet(gc.open("1000GenomesRepeatValidationDatabase").worksheet("RefProductSizes"))

rmloci = ["AR","NOS1","FMR1","CNBP","chr12_131901040_T","chr12_75962280_T","chr12_92269453_T", \
    "chr15_72443969_T","chr16_58599051_A","chr16_67644335_T","chr3_54501407_T", \
    "chr4_176019003_T", "AFF2","ARX","chr7_103989357_CCG","RUNX2","chr12_4096182_TTCC","PHOX2B","HOXD13","PAPBN1"]
ref_prod_sizes = ref_prod_sizes[~ref_prod_sizes["PrimerID"].isin(rmloci)]

############### Write matrix of product sizes ###########

psize_data = {} # PrimerID -> SampleID -> ProductSize
primer_ids = sorted(list(set(ref_prod_sizes["PrimerID"])))
primer_ids = [item for item in primer_ids if item not in rmloci]
samples = list(set(samp_prod_sizes["SampleID"]))

for pr in primer_ids:
    psize_data[pr] = {}
    psize_data[pr]["ref"] = ref_prod_sizes[ref_prod_sizes["PrimerID"]==pr]["ProductSize"].values[0]
    for sm in samples: psize_data[pr][sm] = "NA"
    
for i in range(samp_prod_sizes.shape[0]):
    pr = samp_prod_sizes["PrimerID"].values[i]
    if pr not in primer_ids: continue
    sm = samp_prod_sizes["SampleID"].values[i]
    prod = samp_prod_sizes["ProductSize"].values[i]
    psize_data[pr][sm] = prod

outf = open("TableS2-ProductSizes.csv", "w")
outf.write(",".join(["PrimerID","RefProductSize"] + samples)+"\n")
for pr in primer_ids:
    items = [pr] + [psize_data[pr]["ref"]] + [psize_data[pr][sm] for sm in samples]
    outf.write(",".join([str(item) for item in items])+"\n")
outf.close()
