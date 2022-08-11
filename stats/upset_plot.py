from upsetplot import from_memberships
from upsetplot import plot
from matplotlib import pyplot


methods = [0,0,0,0]
gang_hip = 0
adv_hip = 0
eh_hip = 0
gang_adv = 0
gang_eh = 0
eh_adv = 0
gang_hip_adv = 0
gang_hip_eh = 0
hip_adv_eh = 0
gang_adv_eh = 0
gang_adv_eh_hip = 0

dict = {"1|1|1|1":0, "1|1|1|0":0, "1|1|0|1":0, "1|0|1|1":0, "0|1|1|1":0,
        "1|1|0|0":0, "0|0|1|1":0, "1|0|1|0":0, "0|1|1|0":0, "0|1|0|1":0, "1|0|0|1":0, "0|0|0|1":0, "1|0|0|0":0,"0|1|0|0":0,"0|0|1|0":0 }

for i in range(1,23):
    print(i)
    with open("methods_chr" + str(i) + ".txt") as f:
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
pyplot.savefig("upset_cnt.png", dpi = 1200)

plot(example, show_counts=False, sort_by='cardinality', facecolor="darkblue", shading_color="lightgray")
pyplot.savefig("upset.png", dpi = 1200)



