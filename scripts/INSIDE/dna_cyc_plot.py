import os
import numpy as np
import matplotlib.pyplot as plt

dir40 = "/Users/blake/Downloads/cyc"
fs = os.listdir(dir40)
gens = {"G1": [], "G30": [], "G60": [], "G90": []}
for f in fs:
    f_cell = f.split("_")
    if f_cell[3] in gens:
        f_vals = []
        with open(f"{dir40}/{f}") as file:
            for line in file:
                cell = line.split(",")
                if cell[1] != "C0S_norm":
                    f_vals.append(float(cell[1]))
        gens[f_cell[3]].append(f_vals)
npg0 = np.mean(np.array(gens["G1"]), axis=0)
for g in gens:
    g_int = int(g.split("G")[-1])

    npg = np.mean(np.array(gens[g]), axis=0)
    print(npg.shape)
    plt.plot(npg - npg0, label=g, alpha=g_int / 160 + 0.2, color="black")


plt.legend(title="Generation")

plt.show()
# with open()
