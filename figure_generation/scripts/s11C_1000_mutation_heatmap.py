import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import INSIDE_1000_enhs

# lines.append(f"{chrom1}\t{begin}\t{ending}\n")


def find_differences(str1, str2):
    if len(str1) != len(str2):
        raise ValueError("Strings must have the same length")

    differences = []
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            differences.append(i)

    return differences


def find_consensus_diffs(diff_holder, str1):
    differences = []
    for i, base in enumerate(str1):
        if base not in diff_holder[i]:
            diff_holder[i].append(base)
            differences.append(i)
    return diff_holder, differences


np_heat = np.zeros([300, 1024])


def set_style():
    fig = plt.figure()
    ax = fig.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot
    ax.spines[["right", "top"]].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.5)  # Bottom (x-axis) spine
    ax.spines["left"].set_linewidth(1.5)  # Left (y-axis) spine
    ax.tick_params(axis="both", width=2)


# set_style()

enhs = INSIDE_1000_enhs.return_enhs()
enh_dict = {}
enh_len_lookup = []
for e in enhs:
    chrom1, coords1 = e[1].split(":")
    b1, e1 = coords1.split("-")
    chrom2, coords2 = e[2].split(":")
    b2, e2 = coords2.split("-")
    e1_len = int(e1) - int(b1)
    e2_len = int(e2) - int(b2)
    enh_len_lookup.append((e1_len, e2_len))
enh_strings = np.zeros([1000, 50, 2]).tolist()
with open("INSIDE_1000.txt") as file:
    for line in file:
        cell = line.split("\t")
        enh_num = int(cell[1].split("_")[1])

        gen_num = int(cell[0].split(": ")[-1]) - 1
        e1_len, e2_len = enh_len_lookup[enh_num]
        e1 = cell[2][:e1_len]
        e2 = cell[2][e1_len:]

        enh_strings[enh_num][gen_num][0] = e1
        enh_strings[enh_num][gen_num][1] = e2
        # print(cell[0], len(e1), len(e2))
        # print(cell)
        # print(cell)
    # print(enhs)

# heatmap = []
# cale_size = 180
for scale_size, xticks in [(180, (90, 181, 90)), (10, (5, 11, 5))]:
    current_diff_across = np.zeros([50, scale_size])
    for ei, e in enumerate(enh_strings):
        e1_len, e2_len = enh_len_lookup[ei]

        diff1_holder = []
        for seq_pos in range(e1_len):
            diff1_holder.append([])
        diff2_holder = []
        for seq_pos in range(e2_len):
            diff2_holder.append([])
        for gen in range(0, 50):
            diff1_holder, snp_poses1 = find_consensus_diffs(diff1_holder, e[gen][0])
            diff2_holder, snp_poses2 = find_consensus_diffs(diff2_holder, e[gen][1])
            for p in snp_poses1:
                pos = int(scale_size * (p / e1_len))
                current_diff_across[gen, pos] += 1
            for p in snp_poses2:
                pos = int(scale_size * (p / e2_len))
                current_diff_across[gen, pos] += 1

    plt.figure()
    row_sums = current_diff_across.sum(axis=1, keepdims=True)
    current_diff_across_normalized = current_diff_across / row_sums
    sns.heatmap(current_diff_across_normalized[1:], cmap="Greys")
    # result = np.sum(np_heat, axis=0)
    # print(result.shape)
    # plt.hist(result, bins=128)
    x_tick_positions = np.arange(xticks[0], xticks[1], xticks[2])  # Adjust interval as needed
    plt.xticks(ticks=x_tick_positions, labels=x_tick_positions, rotation=0, ha="center")
    y_tick_positions = np.arange(10, 50, 10)  # Adjust interval as needed
    plt.yticks(ticks=y_tick_positions, labels=y_tick_positions, rotation=0)
    # plt.xlabel("Position (bp)")
    plt.ylabel("Generation")
    # plt.savefig("map_snps.png", dpi=300)
    plt.show()
