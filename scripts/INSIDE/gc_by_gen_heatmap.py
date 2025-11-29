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

gens_to_meta = ["G1", "G50", "G200", "G400"]
heatmaps = {"INSIDE": [], "Shift": []}
for label in ["INSIDE"] + list(range(1, 101)):
    meta_holder = {}
    heatmap = []
    print(label)
    for gen in range(1, 300):
        gen = f"G{gen}"
        left_gen_gc_by_pos = np.zeros(512)
        left_gen_base_by_pos = np.zeros(512)
        right_gen_gc_by_pos = np.zeros(512)
        right_gen_base_by_pos = np.zeros(512)
        for sample in line_seq_holder:
            seq = line_seq_holder[sample][gen][0][0]
            line_sample = int(sample.split("rep")[-1].split("_")[0])
            if label != "INSIDE":
                seq = random_shift(seq, seed=line_sample + label)
            for i, base in enumerate(seq):
                if i >= 512:
                    if base == "G" or base == "C":
                        right_gen_gc_by_pos[i - 512] += 1
                    right_gen_base_by_pos[i - 512] += 1
                else:
                    if base == "G" or base == "C":
                        left_gen_gc_by_pos[i] += 1
                    left_gen_base_by_pos[i] += 1

        window = 1
        gen_gc_by_pos = left_gen_gc_by_pos + right_gen_gc_by_pos
        gen_base_by_pos = right_gen_base_by_pos + left_gen_base_by_pos
        l_gc_ma = np.convolve(left_gen_gc_by_pos / left_gen_base_by_pos, np.ones(window) / window, mode="valid")
        r_gc_ma = np.convolve(right_gen_gc_by_pos / right_gen_base_by_pos, np.ones(window) / window, mode="valid")
        gc_ma = np.convolve(gen_gc_by_pos / gen_base_by_pos, np.ones(window) / window, mode="valid")
        heatmap.append(gc_ma)

        # sns.lineplot(gc_ma, label=gen, alpha=0.8, linewidth=1.25)
    """
    ax = sns.heatmap(data=np.array(heatmap))
    ax.set_xticks([127, 255, 383])  # Python uses 0-based indexing
    ax.set_xticklabels(["128", "256", "384"])
    plt.xticks(rotation=0)
    y_ticks = [1] + list(range(30, 300, 30))
    ax.set_yticks(y_ticks)
    ax.set_yticklabels([str(y) for y in y_ticks])
    plt.xlabel("Position (bp)")
    plt.ylabel("Generation")
    plt.savefig("s8_gc.png", dpi=300)
    plt.show()
    """
    if label == "INSIDE":
        heatmaps[label] = np.array(heatmap)
    else:
        heatmaps["Shift"].append(np.array(heatmap))

heatmaps["Shift"] = np.hstack(heatmaps["Shift"])
for start in range(0, 300, 30):
    plt.figure()
    plt.ylim(0, 1.1)
    plt.xlim(0.3, 1)
    print(label)
    sns.ecdfplot(data=(heatmaps["INSIDE"])[start : start + 30, :].flatten(), label="INSIDE")
    sns.ecdfplot(data=(heatmaps["Shift"])[start : start + 30, :].flatten(), label="Shift")
    plt.title(f"{start}-{start+30}")
    plt.legend()
    plt.savefig(f"gc_ecdf{start}.png")
    plt.show()

plt.axvline(x=128, color="b", linestyle="--")
plt.axvline(x=512 - 128, color="b", linestyle="--")
plt.axhline(y=0.5, color="black", linestyle="--")
plt.title("Moving Average GC Content")
plt.show()
left = []
right = []


# RED Note, shift sequence and wrap back around for each of the 40 lines to see if we get a scramble pattern different from the INSIDE selection

# print(line_seq_holder)
