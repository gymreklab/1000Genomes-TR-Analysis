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

## Adjusting positions for each caller
gang_dict = {
    27573485:27573524
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
    6936729:6936717,
    87604288:87604283,
    63912686:63912685,
    27573484:27573485,
    2442377:2442376,
    1273636:1273621,
    75404064:75404060,
    70011632:70011630,
    7573485:27573529,
    36061394:36061393,
    50595021:50595020,
    18349308:18349307,
    89959098:89959097,
    46950717:46950716,
    76108862:76108861,
    242233155:242233154,
    8128309:8128308,
    80839290:80839285,
    102648430:102648429,
    85908552:85908538,
    92071011:92071009,
    111598951:111598950,
    16327636:16327634,
    103989357:103989356
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
                            gbs = [(len(al) - len(alleles[0]))/period for al in alls]
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
                            gbs = [int(gb)/period for gb in gbs]
                            CNs.append(str(gbs[0]) + "," + str(gbs[1]))
            df.loc[len(df)] = [chrom,pos] + CNs
    df = df.melt(id_vars=['chr','pos'], var_name="Sample", value_name=caller)
    df['pos'] = pd.to_numeric(df['pos'])
    df = pd.merge(df, loci, left_on = "pos", right_on = "Start", how="left")
    df = df[['LocusID','Sample',caller]]
    df.columns = ['PrimerID','SampleID',caller]
    return df.drop_duplicates(keep="first")

# Y = X + offset
def find_offset (X, Y, period):
    useind = [i for i in range(len(X)) if not np.isnan(Y[i]) and not np.isnan(X[i])]
    newX = [X[i] for i in useind]
    newY = [Y[i] for i in useind]

    if len(newX) == 0: return np.nan, np.nan, np.nan
    # measure: diffs in sizes
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
        
        if diff < best_diffs:
            best_diffs = diff
            best_offset = offset
            most_matches = num_matches
    return best_offset, most_matches, best_diffs

def learn_offsets(cap, hipstr_calls, gangstr_calls, loci):
    """
    return df with:
    PrimerId, batch, offset, offset_hipstr, offset_gangstr
    """
    allsamples = list(set(cap["sample"]))

    # Initialize df items
    loci_ = []
    batch_ = []
    offset_ = []
    offset_hipstr_ = []
    offset_gangstr_ = []

    for PrimerID in sorted(list(set(cap["PrimerID"]))):
        period = len(loci[loci["LocusID"]==PrimerID]["Motif"].values[0])
        X = [] # product size from capillary - ref
        Y_gstr = [] # product size from gangstr
        Y_hstr = [] # product size from hipstr

        ref_prod_size = cap[cap["PrimerID"]==PrimerID]["RefProductSize"].values[0]

        for sample in allsamples:
            cap_prod_sizes = cap[(cap["PrimerID"]==PrimerID) & (cap["sample"]==sample)]["Cap"].values[0].split("/")
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

    # Return the df
    return pd.DataFrame({
        "PrimerID": loci_,
        "batch": batch_,
        "offset": offset_,
        "offset_hipstr": offset_hipstr_,
        "offset_gangstr": offset_gangstr_
        })
