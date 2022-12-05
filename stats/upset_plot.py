from upsetplot import from_memberships
from upsetplot import plot
from matplotlib import pyplot
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

dict = {"1|1|1|1":0, "1|1|1|0":0, "1|1|0|1":0, "1|0|1|1":0, "0|1|1|1":0,
        "1|1|0|0":0, "0|0|1|1":0, "1|0|1|0":0, "0|1|1|0":0, "0|1|0|1":0, "1|0|0|1":0, "0|0|0|1":0, "1|0|0|0":0,"0|1|0|0":0,"0|0|1|0":0 }

dir_addr="/projects/ps-gymreklab/helia/ensembl/experiments/upset_plot"
for i in range(1,23):
    print(i)
    with open(dir_addr + "/methods_chr" + str(i) + ".txt") as f:
        for line in f:
            if "Number" in line:
                continue
            else:
                line = line.replace("METHODS=","")
                line = line.replace("\n", "")
                dict[line] += 1

example = from_memberships(
    [['GangSTR'],
    ['HipSTR'],
    ['adVNTR'],
    ['ExpansionHunter'],
    ['HipSTR', 'GangSTR'],
    ['HipSTR', 'adVNTR'],
    ['HipSTR', 'ExpansionHunter'],
    ['GangSTR', 'adVNTR'],
    ['GangSTR', 'ExpansionHunter'],
    ['adVNTR', 'ExpansionHunter'],
    ['adVNTR', 'ExpansionHunter','GangSTR'],
    ['adVNTR', 'ExpansionHunter', 'HipSTR'],
    ['GangSTR', 'ExpansionHunter', 'HipSTR'],
    ['GangSTR', 'adVNTR', 'HipSTR'],
    ['GangSTR', 'adVNTR', 'HipSTR','ExpansionHunter'],
    ],
    data=[dict["0|0|0|1"], dict['0|0|1|0'], dict['1|0|0|0'], dict['0|1|0|0'], dict["0|0|1|1"], dict['1|0|1|0'],
          dict['0|1|1|0'], dict['1|0|0|1'], dict['0|1|0|1'], dict['1|1|0|0'], dict["1|1|0|1"], dict["1|1|1|0"],
          dict["0|1|1|1"], dict["1|0|1|1"] , dict['1|1|1|1']]
    )

plot(example, show_counts="%d", sort_by='cardinality', facecolor="darkblue", shading_color="lightgray")
pyplot.savefig(dir_addr+"/upset_cnt.pdf", dpi = 1200)

plot(example, show_counts=False, sort_by='cardinality', facecolor="darkblue", shading_color="lightgray")
pyplot.savefig(dir_addr+"/upset.pdf", dpi = 1200)



