import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import os
import numpy as np

directory = "../data/mar6_msp/"
# directory = "../data/beds_sep4/"
entries = os.listdir(directory)
"#9cd7a4"
"#c03830"
"#7192BE"


def intersect(region1, region2):
    start1, end1 = region1
    start2, end2 = region2

    start = max(start1, start2)
    end = min(end1, end2)

    if start <= end:
        return (start, end)
    else:
        return None


def set_style(ax):
    # fig = plt.figure(figsize=figsize)
    # ax = fig.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot
    ax.spines[["right", "top", "bottom"]].set_visible(False)
    # ax.spines["bottom"].set_linewidth(1.5)  # Bottom (x-axis) spine
    ax.spines["left"].set_linewidth(1.5)  # Left (y-axis) spine
    # ax.tick_params(axis="both", width=2)
    ax.tick_params(axis="x", width=0, labelbottom=False)
    ax.tick_params(axis="y", width=2, labelsize=6)


def odds_ratio(f):
    a = 0
    b = 0
    c = 0
    d = 0
    size_threshold = 150
    if ".bedgraph" not in f:
        # f = f"{directory}{f}"
        fiber5000 = np.ones([5000])
        coverage5000 = np.zeros([5000])
        with open(f) as file:
            for fi, line in enumerate(file):
                cell = line.split("\t")
                sizes = cell[10]
                starts = cell[11]
                for i in range(int(cell[1]), int(cell[2])):
                    if i < 5000:
                        fiber5000[i - 1] += 1
                p1intersect = False
                p2intersect = False
                for size, start in zip(sizes.split(","), starts.split(",")):
                    size = int(size)
                    start = int(start)

                    if size > size_threshold:
                        block_start = start + int(cell[1])
                        block_end = block_start + size
                        p1 = intersect((2448, 2964), (block_start, block_end))
                        p2 = intersect((2964, 2964 + (512)), (block_start, block_end))
                        # p1 = intersect((1000, 1319), (block_start, block_end))
                        # p2 = intersect((1402, 1872), (block_start, block_end))
                        if p1 != None:
                            p1intersect = True
                        if p2 != None:
                            p2intersect = True
                # print(fi, p1intersect, p2intersect)
                if p1intersect == True:
                    if p2intersect == True:
                        a += 1
                    else:
                        b += 1
                else:
                    if p2intersect == True:
                        c += 1
                    else:
                        d += 1
        try:
            stat1 = round(a**2 / ((b + c) / 2) ** 2, 2)
            stat2 = round((a * d) / (b * c), 2)
            print(f, a, b, c, d, stat1)
        except Exception as e:
            stat1 = -1
            stat2 = -1
        return stat1
        # print(f, a, b, c, d, stat1)


fi_total = 0
filter_entries = []
for f in entries:
    msps = []
    print(f)
    if ".bedgraph" not in f and (".2." in f):
        filter_entries.append(f)
        fi_total += 1
print(fi_total, filter_entries)

filter_entries = [
    "sepr6g60.bam.msp.bed",
    "sep160.bam.msp.bed",
    "sep601.bam.msp.bed",
    "sepr6g1.bam.msp.bed",
    # "sepr6g1_0.bam.msp.bed",
    # ,
]
fe60 = ["marg60.bam.msp.bed", "mar160.bam.msp.bed", "mar601.bam.msp.bed", "marg1.bam.msp.bed"]
fe90 = ["marg90.bam.msp.bed", "mar190.bam.msp.bed", "mar901.bam.msp.bed", "marg1.bam.msp.bed"]

for filter_entries, file_save in ((fe60, "../figures/4c_switch60.svg"), (fe90, "../figures/4c_switch90.svg")):
    """
    filter_entries = [
        "sepr6g90.bam.msp.bed",
        "sep190.bam.msp.bed",
        "sep901.bam.msp.bed",
        "sepr6g1.bam.msp.bed",
    ]



    filter_entries = [
        "G1.2.bam.msp.bed",
        "G30.2.bam.msp.bed",
        "G60.2.bam.msp.bed",
        "G90.2.bam.msp.bed",
        "G120.2.bam.msp.bed",
        "G150.2.bam.msp.bed",
    ]
    """
    # entires =
    f_n = len(filter_entries)
    fig, ax = plt.subplots(2 * f_n, 1, gridspec_kw={"height_ratios": [2, 1] * f_n}, figsize=(4, f_n), sharex=True)
    for fei, f in enumerate(filter_entries):
        msps = []
        # print(f)
        # if ".bedgraph" not in f and "90" in f or "r6g1" in f:
        f_name = f.split(".")[0]
        f = f"{directory}{f}"
        set_style(ax[2 * fei])
        set_style(ax[2 * fei + 1])
        # get coverage
        fiber5000 = np.ones([5000])
        coverage5000 = np.zeros([5000])

        with open(f) as file:
            # print(f)
            for line in file:
                cell = line.split("\t")
                sizes = cell[10]
                starts = cell[11]
                for i in range(int(cell[1]), int(cell[2])):
                    if i < 5000:
                        fiber5000[i - 1] += 1
                for size, start in zip(sizes.split(","), starts.split(",")):
                    size = int(size)
                    start = int(start)

                    if size > 150:
                        block_start = start + int(cell[1])
                        block_end = block_start + size
                        # print(cell[1], cell[2], block_start, block_end)
                        for i in range(block_start, block_end):
                            if i < 5000:
                                coverage5000[i - 1] += 1
        coverage5000 = coverage5000 / fiber5000
        ax[2 * fei].set_ylim(0, 0.9)
        ax[2 * fei].set_yticks([0.9])
        alpha = 0.85
        if fei == 3:
            alpha = 0.35
        stat = odds_ratio(f)
        ax[2 * fei].text(1, 0, f"{stat}", transform=ax[2 * fei].transAxes, ha="center", va="center", fontsize=8)
        ax[2 * fei].bar(range(5000), coverage5000, width=1.0, color="#c03830", alpha=alpha)
        # ax[2 * fei].set_ylabel(f_name, rotation=0)
        ax[2 * fei].set_ylabel(f_name, rotation=0, fontsize="12")
        ax[2 * fei].yaxis.set_label_coords(-0.33, -0.5)
        rect = Rectangle(
            (2448 + 128, 0),
            256,
            1000,
            edgecolor="grey",
            facecolor="grey",
            alpha=0.15,
            linewidth=0.05,
        )
        ax[2 * fei].add_patch(rect)
        rect = Rectangle(
            (2448 + 512 + 128, 0),
            256,
            1000,
            edgecolor="grey",
            facecolor="grey",
            alpha=0.15,
            linewidth=0.05,
        )
        ax[2 * fei].add_patch(rect)
        # get coverage
        fibers = []
        with open(f) as file:
            fi = 0
            for qi, line in enumerate(file):
                cell = line.split("\t")
                sizes = cell[10].split(",")
                starts = cell[11].split(",")
                fiber_len = int(cell[2]) - int(cell[1])
                if fiber_len < 500:
                    continue
                fi += 1
                msp_amt = 0
                for size in sizes:
                    msp_amt += int(size)
                fibers.append([msp_amt, (int(cell[1]), fiber_len, sizes, starts)])
            fibers = fibers[:100]
            fibers.sort(key=lambda x: x[0], reverse=False)
            for fi, fiber in enumerate(fibers):
                rect = Rectangle(
                    (fiber[1][0], fi),
                    fiber[1][1],
                    1,
                    edgecolor="#f0f0f0",
                    facecolor="#f0f0f0",
                    alpha=1,
                    linewidth=1.00,
                )
                ax[2 * fei + 1].add_patch(rect)
                for size, start in zip(fiber[1][2], fiber[1][3]):
                    # print(fi, size, start)
                    if int(size) > 150:
                        rect = Rectangle(
                            (fiber[1][0] + int(start), fi),
                            int(size),
                            1,
                            edgecolor="#c03830",
                            facecolor="#c03830",
                            alpha=1,
                            linewidth=1.5,
                        )
                        ax[2 * fei + 1].add_patch(rect)
                        # regions = ["chr1\t2448\t2964\nchr1\t2964\t3481"]
                        # regions = ["chr1\t2448\t2964\nchr1\t2964\t3481"]
        ax[2 * fei + 1].spines["left"].set_linewidth(0)
        ax[2 * fei + 1].set_xlim(2448 - 103 * 2.5, 3481 + 102 * 3)

        # Set y-axis ticks to just the maximum
        ax[2 * fei + 1].set_yticks([])
        ax[2 * fei + 1].set_ylim(0, fi)
        # fig.tight_layout()
    plt.subplots_adjust(hspace=0.2)
    plt.savefig(file_save)
    # plt.show()
