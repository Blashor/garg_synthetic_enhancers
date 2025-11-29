bed_lines = []
gen_add = 0
survivors = 4
line_seq_holder = {}
import math
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import random


def random_shift(seq, seed=None):
    if seed is not None:
        random.seed(seed)
    length = len(seq)
    split_point = random.randint(0, length - 1)
    return seq[split_point:] + seq[:split_point]


# INSIDE_256_g117.out
with open("INSIDE_jun17.txt") as file:
    current_stat = 0
    gen_tracker = {}
    replicates_mean = {}
    itr = -1
    first_lines = True
    for line in file:
        line = line.strip()
        if "NEW_FILE" in line:
            gen_add += 1

        if "" not in line:
            continue

        if "Gen:" not in line:
            continue
        itr += 1
        surv_num = (itr % survivors) + 1
        cell = line.split("\t")
        cell[0] = "G" + str(int(cell[0].split("Gen:")[-1]) + gen_add)
        current_stat += float(cell[3])
        cell[1] += "_" + str(surv_num)
        if cell[1] not in line_seq_holder:
            line_seq_holder[cell[1]] = {}
        if cell[0] not in line_seq_holder[cell[1]]:
            line_seq_holder[cell[1]][cell[0]] = []
        line_seq_holder[cell[1]][cell[0]].append((cell[2], cell[3]))

gens_to_meta = ["G1", "G90", "G150", "G300"]

meta_holder = {}

for gen in gens_to_meta:
    left_gen_gc_by_pos = np.zeros(512)
    left_gen_base_by_pos = np.zeros(512)
    right_gen_gc_by_pos = np.zeros(512)
    right_gen_base_by_pos = np.zeros(512)
    for sample in line_seq_holder:
        seq = line_seq_holder[sample][gen][0][0]
        line_sample = int(sample.split("rep")[-1].split("_")[0])
        # seq = random_shift(seq, seed=line_sample + 1)
        for i, base in enumerate(seq):
            if i >= 512:
                if base == "G" or base == "C":
                    right_gen_gc_by_pos[i - 512] += 1
                right_gen_base_by_pos[i - 512] += 1
            else:
                if base == "G" or base == "C":
                    left_gen_gc_by_pos[i] += 1
                left_gen_base_by_pos[i] += 1

    window = 10
    gen_gc_by_pos = left_gen_gc_by_pos + right_gen_gc_by_pos
    gen_base_by_pos = right_gen_base_by_pos + left_gen_base_by_pos
    l_gc_ma = np.convolve(left_gen_gc_by_pos / left_gen_base_by_pos, np.ones(window) / window, mode="valid")
    r_gc_ma = np.convolve(right_gen_gc_by_pos / right_gen_base_by_pos, np.ones(window) / window, mode="valid")
    gc_ma = np.convolve(gen_gc_by_pos / gen_base_by_pos, np.ones(window) / window, mode="valid")
    sns.lineplot(gc_ma, label=gen, alpha=0.8, linewidth=1.25)
plt.ylim(0, 1)
plt.axvline(x=128, color="b", linestyle="--")
plt.axvline(x=512 - 128, color="b", linestyle="--")
plt.axhline(y=0.5, color="black", linestyle="--")
plt.title("Moving Average GC Content")
plt.show()
left = []
right = []
# print(line_seq_holder)
