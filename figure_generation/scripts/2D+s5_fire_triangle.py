import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math

# Create a 10x10 matrix filled with zeros


def short_to_small(short_file, locus="chr3:34747631-34762521"):
    chrom, coord_range = locus.split(":")
    start, end = map(int, coord_range.split("-"))
    se_size = end - start
    print(se_size)
    short_lines = []
    needed_size = int((15000 - (se_size)) / 2)
    with open(short_file) as file:
        for line in file:
            cell = line.split("\t")
            if cell[1] == chrom:
                cell[2] = int(cell[2])
                if cell[2] > start and cell[2] < end:
                    short_lines.append(line.strip() + "\n")
    with open(short_file + ".small", "w") as file:
        file.writelines(short_lines)


def short_to_fire_triangle(short_file, locus="chr3:34747631-34762521", resolution=128):
    chrom, coord_range = locus.split(":")
    start, end = map(int, coord_range.split("-"))
    se_size = end - start

    print(se_size)
    bins = list(range(start, start + 15000, 128))
    print(bins)
    matrix = np.zeros((500, 500))
    mat_dict = {}
    with open(short_file) as file:
        for line in file:
            cell = line.split("\t")
            if cell[1] == chrom:
                cell[2] = int(cell[2])
                cell[6] = int(cell[6])
                bin1 = math.floor((cell[2] - start) / resolution)
                bin2 = math.floor((cell[6] - start) / resolution)
                if float(cell[8]) > 10:
                    print(bin1, bin2, cell[8])
                # print(cell[8])
                bins = f"{bin1}_{bin2}"
                if bins not in mat_dict:
                    mat_dict[bins] = []
                mat_dict[bins].append(float(cell[8]))
    for bins in mat_dict:
        bin1, bin2 = bins.split("_")
        matrix[int(bin1), int(bin2)] = np.mean(mat_dict[bins])
    # print(np.log2(matrix[:150, :150]))
    """
    sns.heatmap(
        np.log2(matrix[: int(se_size / resolution), : int(se_size / resolution)]),
        cmap="RdBu_r",
        annot=False,
        cbar=False,
        center=0,
        vmax=6,
        vmin=-6,
        square=True,
    )
    """
    print(start, start + int(se_size / resolution) * resolution, int(se_size / resolution))
    # print(matrix[:150, :150])
    print(matrix[:29, :29].shape)
    sns.heatmap(
        matrix[: int(se_size / resolution), : int(se_size / resolution)],
        cmap=sns.color_palette("Reds", as_cmap=True),
        annot=False,
        cbar=False,
        vmax=5,
        square=True,
    )
    # Show the plot
    plt.title(f"{locus}")
    plt.savefig(f"{locus}.svg")
    plt.show()


"chr2:152396501-152411724"
c_locus = "chr13:48536592-48552015"
c_locus = "chr3:34747631-34762521"
c_locus = "chr3:34745000-34765000"
short_to_small("../data/bin128f_fire_mm10.short", locus=c_locus)
short_to_fire_triangle("../data/bin128f_fire_mm10.short.small", locus=c_locus)
