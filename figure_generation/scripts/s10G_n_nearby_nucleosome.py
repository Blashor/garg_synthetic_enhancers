import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

enhs = []
# lines.append(f"{chrom1}\t{begin}\t{ending}\n")


def in_window(x, period=147, w=2):
    # Fold x into one period
    if x < period - w:
        return False
    r = x % period

    # Check if it's within ±w of the start of the period (0)
    return r <= w or r >= period - w


# for x in [0, 5, 10, 140, 147, 150, 154, 155, 160, 280, 290]:
#    print(x, in_window(x))


from collections import defaultdict


enh_dict = {}

import matplotlib.pyplot as plt
from collections import defaultdict
import bisect

fig, axs = plt.subplots(4, 1, figsize=(6, 8), sharex=True)


seen_coords = set()
for ce_mode in range(2):
    relative_positions = [[], [], [], []]
    relative_fibers = [[], [], [], []]
    true_coords_by_chrom = defaultdict(list)
    enhs = []

    with open("/Users/blake/Documents/gargLab/Figures_for_Box/Figure_1_Ranking/data/3runs_ce_rank.txt") as file:
        for line in file:
            if ce_mode == 0 and "Super" in line:
                cell = line.split("\t")
                if cell[0] in seen_coords:
                    continue  # skip duplicates
                enhs.append((f"{cell[-2]}{cell[-1]}", cell[0], cell[1]))
                chrom, coords = cell[0].split(":")
                chrom2, coords2 = cell[1].split(":")
                tc = (int(coords.split("-")[0]) + int(coords.split("-")[1])) / 2
                # tc = int(coords.split("-")[0])
                true_coords_by_chrom[chrom].append(tc)
                seen_coords.add(cell[0])
            elif ce_mode == 1 and "Super" not in line:
                cell = line.split("\t")
                if cell[0] in seen_coords:
                    continue  # skip duplicates
                enhs.append((f"{cell[-2]}{cell[-1]}", cell[0], cell[1]))
                chrom, coords = cell[0].split(":")
                chrom2, coords2 = cell[1].split(":")
                tc = (int(coords.split("-")[0]) + int(coords.split("-")[1])) / 2
                # tc = int(coords.split("-")[0])
                true_coords_by_chrom[chrom].append(tc)
                seen_coords.add(cell[0])

            else:
                continue

    for chrom in true_coords_by_chrom:
        true_coords_by_chrom[chrom].sort()

    WINDOW = 600
    bins = np.linspace(-WINDOW, WINDOW, int(WINDOW / 2))
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
                block_count = int(cell[9])
                block_sizes = list(map(int, cell[10].rstrip(",").split(",")))
                block_starts = list(map(int, cell[11].rstrip(",").split(",")))

                mut_list = true_coords_by_chrom[chrom]

                # mutations close enough to interact with this nucleosome
                left_c = bisect.bisect_left(mut_list, chrom_start - WINDOW)
                right_c = bisect.bisect_right(mut_list, chrom_end + WINDOW)

                nuc_starts = []
                nuc_ends = []
                for i in range(block_count):
                    size = block_sizes[i]
                    if size <= 140:
                        continue

                    nuc_start = chrom_start + block_starts[i]
                    nuc_end = nuc_start + size
                    nuc_center = (nuc_start + nuc_end) // 2

                    # mutations close enough to interact with this nucleosome
                    left = bisect.bisect_left(mut_list, nuc_start - WINDOW)
                    right = bisect.bisect_right(mut_list, nuc_end + WINDOW)

                    for tc in mut_list[left:right]:
                        # relative positions of this mutation to the nucleosome span

                        start_pos = max(nuc_start, tc - WINDOW)
                        end_pos = min(nuc_end, tc + WINDOW)
                        # print(start_pos, end_pos, start_pos - end_pos)
                        # convert to relative coordinates

                        footprint = np.arange(start_pos, end_pos, dtype=int) - tc
                        f_footprint = np.arange(chrom_start, chrom_end, dtype=int) - tc
                        # relative_positions.append(np.arange(start_pos, end_pos, dtype=int) - tc)
                        nuc_starts.append((nuc_center - tc, footprint, f_footprint))
                        # print((nuc_center - tc, footprint))
                if len(nuc_starts) > 3:
                    greatest_nonpositive = sorted(nuc_starts, key=lambda x: abs(x[0]))
                    # greatest_nonpositive = sorted((t for t in nuc_starts if t[0] <= 0), key=lambda x: abs(x[0]))
                    if len(greatest_nonpositive) > 3:
                        for phase in range(4):
                            # greatest_nonpositive = min(nuc_starts, key=lambda x: abs(x[0]))
                            #
                            if greatest_nonpositive != None:
                                relative_fibers[phase].append(greatest_nonpositive[phase][2])
                                relative_positions[phase].append(greatest_nonpositive[phase][1])
    # after the loop
    for phase in range(4):
        relative_positions[phase] = np.concatenate(relative_positions[phase])
        relative_fibers[phase] = np.concatenate(relative_fibers[phase])

        ax = axs[phase]
        ax.spines[["right", "top"]].set_visible(False)
        ax.spines["bottom"].set_linewidth(1.5)  # Bottom (x-axis) spine
        ax.spines["left"].set_linewidth(1.5)  # Left (y-axis) spine
        ax.tick_params(axis="both", width=2)
        ax.set_ylim([0, 90])
        nuc_hist, _ = np.histogram(relative_positions[phase], bins=bins)
        # ratio = nuc_hist / min(nuc_hist)
        fiber_hist, _ = np.histogram(relative_fibers[phase], bins=bins)

        fiber_hist = np.where(fiber_hist == 0, 1, fiber_hist)
        ratio = 100 * nuc_hist / fiber_hist

        bin_centers = (bins[:-1] + bins[1:]) / 2
        bin_width = bins[1] - bins[0]
        # sns.ecdfplot(relative_positions, label=ce_mode)
        # ax.bar(bin_centers, ratio, width=bin_width, align="center", alpha=0.6)
        ax.plot(bin_centers, ratio)
    # dratio_dx = np.gradient(ratio, bin_width)
    # ax.set_ylim([-0.45, 0.45])
    # ax.plot(bin_centers, dratio_dx, linewidth=2)
    # for i in range(1, 7):
    #    ax.axvline(i * -147, linestyle="--")
    #    ax.axvline(i * 147, linestyle="--")
    # if ce_mode == 0:
    #    ax.set_title(f"Super")
    # if ce_mode == 1:
    #    ax.set_title(f"CE")
    # ax.set_xlabel("Nucleosome Occupancy (bp)")

# axs[0].set_ylabel("Mutation Percent")
plt.legend()
plt.tight_layout()
plt.show()
