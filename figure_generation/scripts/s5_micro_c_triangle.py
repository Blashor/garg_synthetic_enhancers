import numpy as np
import hicstraw
import json
import matplotlib.pyplot as plt
import seaborn as sns
import scipy

locus = "chr2:152396501-152411724"
# locus = "chr13:48536592-48552015"
locus = "chr3:34747631-34762521"
locus = "chr3:34748000-34762000"
locus = "chr3:34745000-34765000"


def main():
    micro_hic = hicstraw.HiCFile("../data/4DNFI6HG4GP3.hic")
    current_chrom = ""

    x = []
    y = []
    s = []

    chrom, coord_range = locus.split(":")
    start, end = map(int, coord_range.split("-"))

    chrom_n = chrom.split("chr")[1]
    mzd_micro = micro_hic.getMatrixZoomData(chrom_n, chrom_n, "observed", "VC_SQRT", "BP", 1000)
    micro_mat = mzd_micro.getRecordsAsMatrix(start, end, start, end)
    sns.heatmap(data=np.array(micro_mat)[1:, 1:], cmap=sns.color_palette("Reds", as_cmap=True), square=True, vmax=50)
    plt.show()


main()
