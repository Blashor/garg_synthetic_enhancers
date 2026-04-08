import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import INSIDE_1000_enhs

# lines.append(f"{chrom1}\t{begin}\t{ending}\n")


def string_diff(a, b):
    diffs = []
    for i in range(len(a)):
        if a[i] != b[i]:
            diffs.append(i)
    return diffs


def set_style():
    fig = plt.figure()
    ax = fig.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot
    ax.spines[["right", "top"]].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.5)  # Bottom (x-axis) spine
    ax.spines["left"].set_linewidth(1.5)  # Left (y-axis) spine
    ax.tick_params(axis="both", width=2)


set_style()

enhs = INSIDE_1000_enhs.return_enhs()
enh_dict = {}
gen_to_count = 10
for gen_to_count in [10, 20, 50]:
    print(gen_to_count)
    with open("../data/INSIDE_1000.txt") as file:
        for line in file:
            cell = line.split("\t")
            if cell[1] not in enh_dict:
                enh_dict[cell[1]] = {0: "", 1: ""}
            if cell[0] == f"Gen: {gen_to_count-9}":
                enh_dict[cell[1]][0] = cell[2]
            if cell[0] == f"Gen: {gen_to_count}":
                enh_dict[cell[1]][1] = cell[2]
    for e in enhs:
        chrom1, coords1 = e[1].split(":")
        b1, e1 = coords1.split("-")
        chrom2, coords2 = e[2].split(":")
        b2, e2 = coords2.split("-")
        coords = [chrom1, b1, e1, b2, e2]
        enh_dict[f"{e[0]}_rep1"]["coords"] = coords

    from collections import defaultdict

    true_coords_by_chrom = defaultdict(list)

    for enh in enh_dict:
        diffs = string_diff(enh_dict[enh][0], enh_dict[enh][1])

        chrom, b1, e1, b2, e2 = enh_dict[enh]["coords"]

        b1, e1, b2, e2 = map(int, (b1, e1, b2, e2))
        total_len = e2 - b1
        # diffs = list(range(0, total_len, 5))
        # diffs = [(b1 + e1) / 2 - b1, (b2 + e2) / 2 - b1]
        enh1_len = e1 - b1

        for d in diffs:
            if d < enh1_len:
                tc = b1 + d
            else:
                tc = b2 + (d - enh1_len)

            true_coords_by_chrom[chrom].append(tc)
    for chrom in true_coords_by_chrom:
        true_coords_by_chrom[chrom].sort()

    import bisect

    WINDOW = 300
    bins = np.linspace(-WINDOW, WINDOW, 150)
    #
    relative_positions = []
    fiber_positions = []
    files = [
        "../data/1000_D01_2_nuc.bed",
        "../data/1000_D01_nuc.bed",
        "../data/1000_A01_nuc.bed",
    ]
    for f in files:
        with open(f) as file:
            for line in file:
                cell = line.rstrip().split("\t")

                chrom = cell[0]

                if chrom not in true_coords_by_chrom:
                    continue

                chrom_start = int(cell[1])
                chrom_end = int(cell[2])
                chrom_center = (chrom_start + chrom_end) // 2
                block_count = int(cell[9])
                block_sizes = list(map(int, cell[10].rstrip(",").split(",")))
                block_starts = list(map(int, cell[11].rstrip(",").split(",")))

                mut_list = true_coords_by_chrom[chrom]
                left = bisect.bisect_left(mut_list, chrom_start - WINDOW)
                right = bisect.bisect_right(mut_list, chrom_end + WINDOW)

                for tc in mut_list[left:right]:
                    pos = tc - chrom_center
                    if -WINDOW <= pos < WINDOW:
                        fiber_footprint = np.arange(chrom_start, chrom_end, dtype=int) - tc
                        fiber_footprint = fiber_footprint[(fiber_footprint >= -WINDOW) & (fiber_footprint < WINDOW)]
                        fiber_positions.append(fiber_footprint)
                for i in range(block_count):
                    nuc_start = chrom_start + block_starts[i]
                    nuc_end = nuc_start + block_sizes[i]
                    nuc_center = (nuc_start + nuc_end) // 2

                    left = bisect.bisect_left(mut_list, nuc_center - WINDOW)
                    right = bisect.bisect_right(mut_list, nuc_center + WINDOW)

                    for tc in mut_list[left:right]:
                        pos = tc - nuc_center
                        if -WINDOW <= pos < WINDOW:
                            # relative_positions.append(tc)
                            footprint = np.arange(nuc_start, nuc_end, dtype=int) - tc
                            # for pos in range(nuc_start, nuc_end):
                            footprint = footprint[(footprint >= -WINDOW) & (footprint < WINDOW)]
                            if len(footprint) > 1:
                                relative_positions.append(footprint)
                            relative_positions.append(footprint)

    relative_positions = np.concatenate(relative_positions)
    fiber_positions = np.concatenate(fiber_positions)
    # counts, edges, patches = plt.hist(relative_positions, bins=150, alpha=0)

    # Then np.histogram uses the same edges
    nuc_hist, _ = np.histogram(relative_positions, bins=bins)
    fiber_hist, _ = np.histogram(fiber_positions, bins=bins)
    # plt.clf()
    fiber_hist = np.where(fiber_hist == 0, 1, fiber_hist)
    # Bar plot matches
    bin_centers = (bins[:-1] + bins[1:]) / 2
    bin_width = bins[1] - bins[0]
    plt.plot(bin_centers, nuc_hist / fiber_hist, alpha=0.8, label=f"{gen_to_count}")
    # sns.ecdfplot(relative_positions)
    plt.ylim([0, 30])
    plt.axvline(-73, linestyle="--")
    plt.axvline(73, linestyle="--")
    plt.xlabel("Distance from nucleosome center (bp)")
    plt.ylabel("Mutation Percent")

plt.legend()
# plt.title("Mutation distribution around nucleosomes")
plt.show()
