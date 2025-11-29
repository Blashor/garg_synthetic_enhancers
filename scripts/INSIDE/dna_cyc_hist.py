import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

mpl.use("Agg")

dir40 = "/home/bmt26/garg/blake/bendable"
fs = os.listdir(dir40)
gens = {0: [], 30: [], 60: [], 90: [], 120: [], 150: []}
gens = {0: [], 10: [], 20: [], 30: [], 40: [], 50: []}
for f in fs:
    if ".txt" in f:
        # print(f)
        f_cell = f.split("_")
        gen_n = int(f_cell[-2].replace("G", ""))
        gen_modulo = (gen_n // 10) * 10 + 10
        # print(gen_n, gen_modulo)
        if gen_modulo not in gens:
            gens[gen_modulo] = []
        with open(f"{dir40}/{f}") as file:
            for line in file:
                cell = line.split(",")
                if cell[1] != "C0S_norm" and cell[1] != "C0_norm":
                    if gen_n == 1:
                        gens[gen_n - 1].append(float(cell[1]))
                    else:
                        gens[gen_modulo].append(float(cell[1]))

    # gens[f_cell[3]].append(f_vals)

for g in gens:
    print(g)
    if g != 180:
        # g_int = g * 30  # int(g.split("G")[-1])
        # npg = np.mean(np.array(gens[g]), axis=0)
        # print(npg.shape)
        if g == 0:
            sns.ecdfplot(np.array(gens[g]), label=g, alpha=1, color="black", linestyle="--")
        else:
            sns.ecdfplot(np.array(gens[g]), label=g, alpha=g / 167 + 0.1)


plt.legend(title="Generation")
plt.xlim([-2, 2])
plt.savefig("/home/bmt26/garg/blake/1000_cyc.png")
# with open()
