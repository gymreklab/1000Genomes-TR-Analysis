import pandas as pd
import numpy as np
import re

pg_samples = [
    "NA12889",
    "NA12890",
    "NA12877",
    "NA12891",
    "NA12892",
    "NA12878",
    "NA12879",
    "NA12880",
    "NA12881",
    "NA12882",
    "NA12883",
    "NA12884",
    "NA12885",
    "NA12886",
    "NA12887",
    "NA12888",
    "NA12893"
]

coriell_samples = [
    "NA24385",
    "NA24149",
    "NA24143",
    "NA24631",
    "NA24694",
    "NA24695",
    "HG00268",
    "NA20503",
    "NA20504",
    "NA20505",
    "NA20506",
    "NA18488",
    "HG02970",
    "NA19434",
    "NA19238",
    "NA19239",
    "NA18939",
    "NA18940",
    "NA18941",
    "NA18942",
    "HG00766",
    "HG03736",
    "NA20846",
    "NA20847",
    "NA20849",
    "NA20862",
    "HG01112",
    "HG01113",
    "HG01119",
    "HG00731",
    "HG00732"
]

## Adjusting positions for each caller
gang_dict = {
    63912686:63912685,
    87604288:87604283,
    6936729:6936717,
    27573524:27573529
}

hip_dict = {
    1273621:1273636,
    75404060:75404064,
    70011630:70011632,
    85908538:85908552,
    92071009:92071011,
    111598950:111598951,
    16327634:16327636,
    27573485:27573529
}

ensemble_dict = {
    1273621:1273636,
    36061393:36061394,
    50595020:50595021,
    89959097:89959098,
    46950716:46950717,
    6936729:6936717,
    76108861:76108862,
    87604288:87604283,
    63912686:63912685,
    27573484:27573485,
    2442376:2442377,
    75404060:75404064,
    70011630:70011632,
    7573485:27573529,
    18349307:18349308,
    242233154:242233155,
    8128308:8128309,
    80839285:80839290,
    102648429:102648430,
    85908538:85908552,
    92071009:92071011,
    111598950:111598951,
    16327634:16327636,
    103989357:103989356,
    27573485:27573529
}

def load_vcf(addr, samples, caller, dict_, loci):
    df = pd.DataFrame(columns=['chr','pos'] + samples)
    with open(addr) as f:
        for line in f:
            line = line.strip().split()
            chrom = line[0]
            pos = int(line[1])
            if pos in dict_:
                pos = dict_[pos]
            period = line[7].split(";")[2]
            period = int(period.replace("PERIOD=",""))
            CNs = []
            alleles = [line[3]] + line[4].split(",")
            if caller == "GangSTR":                
                for call in line[9:]:
                    gts = call.split(":")[0]
                    if gts == ".":
                        CNs.append(".")
                    else:
                        gts = gts.split("/")
                        if pos == 3074877: ## Correcting HTT call
                            ref_cn = alleles[0].count("CAG")
                            gt1_cn = alleles[int(gts[0])].count("CAG")
                            gt2_cn = alleles[int(gts[1])].count("CAG")
                            CNs.append(str(gt1_cn - ref_cn) + "," + str(gt2_cn - ref_cn))
                        else:
                            gts = [int(gt) for gt in gts]
                            alls = [alleles[gt] for gt in gts]
                            gbs = [int((len(al) - len(alleles[0]))/period) for al in alls]
                            CNs.append(str(gbs[0]) + "," + str(gbs[1]))
            else:    
                for call in line[9:]:
                    gbs = call.split(":")[1]
                    gts = call.split(":")[0]
                    if gbs == ".":
                        CNs.append(".")
                    else:
                        if caller == "HipSTR":
                            splitter = "\|"
                        else:
                            splitter = "/"
                        gbs = re.split(splitter,gbs)
                        gts = re.split(splitter,gts)
                        if pos == 3074877: ## Correcting HTT call
                            ref_cn = alleles[0].count("CAG")
                            gt1_cn = alleles[int(gts[0])].count("CAG")
                            gt2_cn = alleles[int(gts[1])].count("CAG")
                            CNs.append(str(gt1_cn - ref_cn) + "," + str(gt2_cn - ref_cn))
                        else:
                            gbs = [int(int(gb)/period) for gb in gbs]
                            CNs.append(str(gbs[0]) + "," + str(gbs[1]))
            df.loc[len(df)] = [chrom,pos] + CNs
    df = df.melt(id_vars=['chr','pos'], var_name="Sample", value_name=caller)
    df['pos'] = pd.to_numeric(df['pos'])
    df = pd.merge(df, loci, left_on = "pos", right_on = "Start", how="left")
    df = df[['LocusID','Sample',caller]]
    df.columns = ['PrimerID','SampleID',caller]
    df[caller] = [OrderAlleles(item) for item in list(df[caller])]
    return df.drop_duplicates(keep="first")

def OrderAlleles(item):
    """
    Ensure the smaller allele is listed first
    """
    alleles = item.split(",")
    if len(alleles) != 2 or "." in item: return "./."
    a1 = int(alleles[0])
    a2 = int(alleles[1])
    if a1 < a2: return "%s,%s"%(a1, a2)
    else: return "%s,%s"%(a2, a1)

# Y = X + offset
def find_offset (X, Y, period):
    useind = [i for i in range(len(X)) if not np.isnan(Y[i]) and not np.isnan(X[i])]
    newX = [X[i] for i in useind]
    newY = [Y[i] for i in useind]

    if len(newX) == 0: return np.nan, np.nan, np.nan
    most_matches = -10
    best_diffs = 100000
    best_offset = -1000
    for offset in range(-30,30):
        x_plus_off = [x + offset for x in newX]
        
        # Get matches after round to nearest repeat unit
        rounded_x = [round(item/period) for item in x_plus_off]
        rounded_y = [round(item/period) for item in newY]
        num_matches = sum(np.equal(rounded_x, rounded_y))
        
        # Get diffs. Take top %percentile since we don't know how many calls are wrong
        diffs = [abs(x_plus_off[i]-newY[i]) for i in range(len(x_plus_off))]
        diff = np.max(sorted(diffs)[0:10]) #np.median(diffs)

#        if diff < best_diffs:
        if num_matches > most_matches:
            best_diffs = diff
            best_offset = offset
            most_matches = num_matches
    return best_offset, most_matches, best_diffs

def learn_offsets(cap, hipstr_calls, gangstr_calls, loci):
    """
    return df with:
    PrimerId, batch, offset, offset_hipstr, offset_gangstr, period
    """
    allsamples = list(set(cap["sample"]))

    # Initialize df items
    loci_ = []
    batch_ = []
    offset_ = []
    period_ = []
    offset_hipstr_ = []
    offset_gangstr_ = []

    for PrimerID in sorted(list(set(cap["PrimerID"]))):
        period = len(loci[loci["LocusID"]==PrimerID]["Motif"].values[0])
        X = [] # product size from capillary - ref
        Y_gstr = [] # product size from gangstr
        Y_hstr = [] # product size from hipstr

        ref_prod_size = cap[cap["PrimerID"]==PrimerID]["RefProductSize"].values[0]

        for sample in allsamples:
            cap_prod_sizes = cap[(cap["PrimerID"]==PrimerID) & (cap["sample"]==sample)]["Prd"].values[0].split("/")
            # note, this filters several SCA2 loci where we saw 3 peaks
            if len(cap_prod_sizes) != 2:
                X.extend([np.nan, np.nan])
            else:
                X.extend(sorted([float(i) - float(ref_prod_size) for i in cap_prod_sizes]))
            
            try:
                gb_gstr = sorted(gangstr_calls[(gangstr_calls["PrimerID"]==PrimerID) & \
                                               (gangstr_calls["SampleID"]==sample)]["GangSTR"].values[0].split(","))
                gb_gstr = [float(item) for item in gb_gstr]
            except:
                gb_gstr = [np.nan, np.nan]
            Y_gstr.extend(gb_gstr)
            
            try:
                gb_hstr = sorted(hipstr_calls[(hipstr_calls["PrimerID"]==PrimerID) & \
                                               (hipstr_calls["SampleID"]==sample)]["HipSTR"].values[0].split(","))
                gb_hstr = [float(item) for item in gb_hstr]
            except:
                gb_hstr = [np.nan, np.nan]
            Y_hstr.extend(gb_hstr)

        offset_gstr, num_matches_gstr, diff_gstr = find_offset(X, Y_gstr, period)
        offset_hstr, num_matches_hstr, diff_hstr = find_offset(X, Y_hstr, period)

        pg_ind = [i for i in range(len(allsamples)) if allsamples[i] in pg_samples]
        X_PG = [X[i] for i in pg_ind]
        Y_gstr_PG = [Y_gstr[i] for i in pg_ind]
        Y_hstr_PG = [Y_hstr[i] for i in pg_ind]
        offset_gstr_PG, num_matches_gstr_PG, diff_gstr_PG = find_offset(X_PG, Y_gstr_PG, period)
        offset_hstr_PG, num_matches_hstr_PG, diff_hstr_PG = find_offset(X_PG, Y_hstr_PG, period)

        coriell_ind = [i for i in range(len(allsamples)) if allsamples[i] not in pg_samples]
        X_Coriell = [X[i] for i in coriell_ind]
        Y_gstr_Coriell = [Y_gstr[i] for i in coriell_ind]
        Y_hstr_Coriell = [Y_hstr[i] for i in coriell_ind]
        offset_gstr_Coriell, num_matches_gstr_Coriell, diff_gstr_Coriell = find_offset(X_Coriell, Y_gstr_Coriell, period)
        offset_hstr_Coriell, num_matches_hstr_Coriell, diff_hstr_Coriell = find_offset(X_Coriell, Y_hstr_Coriell, period)
        
        # Get overall offset
        if not np.isnan(offset_hstr):
            offset = offset_hstr
        elif not np.isnan(offset_gstr):
            offset = offset_gstr
        else: offset = np.nan

        loci_.extend([PrimerID]*2)
        batch_.extend(["PG","Coriell"])
        offset_.extend([offset]*2)
        offset_hipstr_.extend([offset_hstr_PG, offset_hstr_Coriell])
        offset_gangstr_.extend([offset_gstr_PG, offset_gstr_Coriell])
        period_.extend([period]*2)

    # Return the df
    return pd.DataFrame({
        "PrimerID": loci_,
        "batch": batch_,
        "offset": offset_,
        "offset_hipstr": offset_hipstr_,
        "offset_gangstr": offset_gangstr_,
        "period": period_
        })

def GetSimpleBin(minval, maxval, period, start, psize):
    al = start
    for i in range(minval, maxval, period):
        if psize >= i and psize < i+period: return al
        al += 1
    return "."

def GetBinnedAlleles(primer_id, product_sizes, sampleid):
    psizes = [float(item) for item in product_sizes.split("/")]
    alleles = []
    for ps in psizes:
        alleles.append(GetBinnedSingleAllele(primer_id, ps, sampleid))
    if "." in alleles: return "."
    return "%s,%s"%(alleles[0], alleles[1])

def GetBinnedSingleAllele(primer_id, psize, sampleid):
    if primer_id == "ATN1":
        return GetSimpleBin(129, 190, 3, -7, psize)
    if primer_id == "ATXN10":
        return GetSimpleBin(192, 240, 5, -3, psize)
    if primer_id == "CACNA1A":
        return GetSimpleBin(131, 170, 3, -7, psize)
    if primer_id == "chr1_106950468_GT":
        if psize < 197:
            return GetSimpleBin(177, 197, 2, -6, psize)
        else: return GetSimpleBin(197, 220, 2, 3, psize)
    if primer_id == "chr1_217937145_AAAT":
        return GetSimpleBin(174, 200, 4, -5, psize)
    if primer_id == "chr1_242233155_AATG":
        return GetSimpleBin(244, 300, 4, -1, psize)
    if primer_id == "chr1_6733191_CA":
        return GetSimpleBin(168, 200, 2, -1, psize)
    if primer_id == "chr1_76108862_GT":
        return GetSimpleBin(244, 300, 2, -11, psize)
    if primer_id == "chr10_57337770_CA":
        return GetSimpleBin(196, 250, 2, -7, psize)
    if primer_id == "chr10_72427734_AT":
        return GetSimpleBin(198, 250, 2, -10, psize)
    if primer_id == "chr11_2442377_AC":
        return GetSimpleBin(171, 250, 2, -5, psize)
    if primer_id == "chr11_36246191_AATA":
        return GetSimpleBin(191, 250, 4, -3, psize)
    if primer_id == "chr11_63798227_AC":
        return GetSimpleBin(234, 300, 2, -5, psize)
    if primer_id == "chr11_85908552_TAGA":
        return GetSimpleBin(204, 300, 4, -5, psize)
    if primer_id == "chr12_70011632_AG":
        if sampleid in coriell_samples:
            return GetSimpleBin(226, 260, 2, -4, psize)
        else: return GetSimpleBin(225, 260, 2, -4, psize)
    if primer_id == "chr13_102648430_GT":
        return GetSimpleBin(214, 250, 2, -7, psize)
    if primer_id == "chr13_75404064_AC":
        return GetSimpleBin(228, 280, 2, -4, psize)
    if primer_id == "chr13_81527591_CTAT":
        return GetSimpleBin(296, 350, 4, -1, psize)
    if primer_id == "chr15_53480481_AC":
        return GetSimpleBin(299, 350, 2, -7, psize)
    if primer_id == "chr15_80839290_GT":
        return GetSimpleBin(203, 250, 2, 0, psize)
    if primer_id == "chr16_62573638_AATA":
        return GetSimpleBin(249, 300, 4, -2, psize)
    if primer_id == "chr17_8128309_TAGA":
        return GetSimpleBin(380, 420, 4, -3, psize)
    if primer_id == "chr18_38020886_TCTA":
        return GetSimpleBin(380, 450, 4, -5, psize)
    if primer_id == "chr2_1273636_TA":
        return GetSimpleBin(234, 280, 2, -1, psize)
    if primer_id == "chr20_18095896_AC":
        return GetSimpleBin(178, 250, 2, -9, psize)
    if primer_id == "chr20_42508128_TG":
        return GetSimpleBin(233, 280, 2, 0, psize)
    if primer_id == "chr3_46654755_TG":
        return GetSimpleBin(312, 350, 2, -2, psize)
    if primer_id == "chr4_123573491_TTAT":
        return GetSimpleBin(201, 250, 4, -1, psize)
    if primer_id == "chr4_136965932_AT":
        return GetSimpleBin(176, 250, 2, -3, psize)
    if primer_id == "chr4_46950717_AAAT":
        return GetSimpleBin(304, 350, 4, 0, psize)
    if primer_id == "chr5_18768193_GT":
        if psize < 177:
            return GetSimpleBin(157, 177, 2, -8, psize)
        else: return GetSimpleBin(176, 250, 2, 1, psize)
    if primer_id == "chr6_12322154_AC":
        return GetSimpleBin(258, 300, 2, -3, psize)
    if primer_id == "chr6_50207587_TCTA":
        return GetSimpleBin(400, 450, 4, -2, psize)
    if primer_id == "chr6_84874972_AC":
        return GetSimpleBin(239, 300, 2, -2, psize)
    if primer_id == "chr6_89959098_CA":
        return GetSimpleBin(210, 300, 2, -11, psize)
    if primer_id == "chr7_18349308_AC":
        if sampleid in coriell_samples:
            return GetSimpleBin(207, 300, 2, -7, psize)
        else: return GetSimpleBin(206, 300, 2, -7, psize)
    if primer_id == "chr7_27264534_AC":
        return GetSimpleBin(224, 250, 2, -4, psize)
    if primer_id == "chr8_50595021_TG":
        return GetSimpleBin(224, 275, 2, -3, psize)
    if primer_id == "chr9_36061394_AC":
        return GetSimpleBin(250, 300, 2, -6, psize)
    if primer_id == "chr9_6685998_TTA":
        return GetSimpleBin(265, 300, 3, 0, psize)
    if primer_id == "DMPK":
        if sampleid in coriell_samples:
            return GetSimpleBin(855, 950, 3, -17, psize)
        else: return GetSimpleBin(856, 950, 3, -15, psize)
    if primer_id == "HTT":
        return GetSimpleBin(101, 200, 3, -5, psize)
    if primer_id == "JPH3":
        return GetSimpleBin(113, 200, 3, -6, psize)
    if primer_id == "PPP2R2B":
        if sampleid in coriell_samples:
            return GetSimpleBin(159, 200, 3, -2, psize)
        else: return GetSimpleBin(162, 200, 3, 0, psize)
    if primer_id == "SCA1":
        return GetSimpleBin(211, 260, 3, -7, psize)
    if primer_id == "SCA2":
        return GetSimpleBin(88, 250, 3, -11, psize)
    if primer_id == "SCA3":
        return GetSimpleBin(114, 200, 3, 0, psize)
    if primer_id == "SCA7":
        if sampleid in coriell_samples:
            return GetSimpleBin(290, 350, 3, -1, psize)
        else: return GetSimpleBin(291, 350, 3, 0, psize)
    return "."

def GetCap(x):
    if x["PrimerID"]=="FXN": return "./."
    cap_prod_sizes = [float(item)-x["RefProductSize"] for item in x["Prd"].split("/")]
    if len(cap_prod_sizes) != 2: return "./."
    return "%s,%s"%(int(round((cap_prod_sizes[0]+x["offset"])/x["period"])), \
                    int(round((cap_prod_sizes[1]+x["offset"])/x["period"])))

def GetMatch(calls1, calls2):
    calls1 = list(calls1)
    calls2 = list(calls2)
    matches = []
    for i in range(len(calls1)):
        if type(calls1[i]) == float or type(calls2[i]) == float:
            matches.append(np.nan)
        elif "." in calls1[i] or "." in calls2[i]:
            matches.append(np.nan)
        else: matches.append(calls1[i] == calls2[i])
    return matches
