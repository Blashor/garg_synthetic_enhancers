import math
import sys
import json
import scipy
import numpy
from matplotlib import pyplot as plt
import seaborn as sns
import mpl_style
import numpy as np


#
# Median odds ratio of FIRE pairs at given distance bins
#


def numPts_below_line(myVector, slope, x):
    yPt = myVector[x]
    b = yPt - (slope * x)
    xPts = numpy.arange(1, len(myVector) + 1)
    return numpy.sum(myVector <= (xPts * slope + b))


def calculate_cutoff(inputVector, drawPlot=False, **kwargs):
    inputVector = numpy.sort(inputVector)
    inputVector[inputVector < 0] = 0

    slope = (numpy.max(inputVector) - numpy.min(inputVector)) / len(inputVector)

    result = numpy.vectorize(lambda x: numPts_below_line(inputVector, slope, x))(numpy.arange(len(inputVector)))
    # print(inputVector)
    # print(result)
    xPt = numpy.floor(numpy.argmin(result)) + 1
    # print(xPt)
    y_cutoff = inputVector[int(xPt - 1)]

    if drawPlot:
        plt.plot(numpy.arange(1, len(inputVector) + 1), inputVector, "b-", **kwargs)
        b = y_cutoff - (slope * xPt)
        plt.axvline(x=xPt, color="gray", linestyle="dashed")
        plt.axhline(y=y_cutoff, color="gray", linestyle="dashed")
        plt.scatter(xPt, y_cutoff, marker="o", color="red")
        plt.plot(
            numpy.arange(1, len(inputVector) + 1),
            slope * numpy.arange(1, len(inputVector) + 1) + b,
            "r-",
        )
        plt.title(
            f"x={xPt}\ny={y_cutoff:.3f}\nFold over Median={y_cutoff / numpy.median(inputVector):.3f}x\nFold over Mean={y_cutoff / numpy.mean(inputVector):.3f}x"
        )
        plt.axvline(x=numpy.sum(inputVector == 0), color="pink", linestyle="dotted", linewidth=2)
        plt.show()
    return y_cutoff


def prep_object():
    lines = []
    with open(prefix + "_Cov.bed", "r") as file:
        lines = file.read().split("\n")
    for line in lines:
        fiber_id = line[: line.rfind("\t")]
        fiber_to_se(line)

    # This moves it into a format that gets total local fiber count so a data matrix can be made
    matrix = []
    for se_id, se in se_fire_score.items():
        if se_id != "":
            se_mat = {"seId": se_id, "enhs": []}
            for enh_id, enh in se.items():
                if enh_id != "fiberList":
                    enh_column = {"enhId": enh_id, "fibers": []}
                    for fiber in se["fiberList"]:
                        try:
                            fire_score = enh[fiber]
                            enh_column["fibers"].append(float(fire_score))
                        except KeyError:
                            enh_column["fibers"].append(None)
                    se_mat["enhs"].append(enh_column)
            matrix.append(se_mat)
    with open(prefix + "_obj.json", "w") as file:
        file.write(json.dumps(matrix))
    return matrix


def fiber_to_se(fiber):
    cells = fiber.split("\t")
    if len(cells) == 8:
        enh_id = cells[0] + ":" + cells[1] + "-" + cells[2]
        super_id = cells[3] + ":" + cells[4] + "-" + cells[5]
        fiber_id = cells[6]
        fire_score = cells[7]
        # print(cells)
        if super_id not in se_fire_score:
            se_fire_score[super_id] = {"fiberList": []}
        if enh_id not in se_fire_score[super_id]:
            se_fire_score[super_id][enh_id] = {}
        if fiber_id not in se_fire_score[super_id]["fiberList"]:
            se_fire_score[super_id]["fiberList"].append(fiber_id)
        # this is to deal with nucleosomes + FIRE elements both overlapping underneath a peak
        if fiber_id in se_fire_score[super_id][enh_id]:
            if fire_score < se_fire_score[super_id][enh_id][fiber_id]:
                se_fire_score[super_id][enh_id][fiber_id] = fire_score
        else:
            se_fire_score[super_id][enh_id][fiber_id] = fire_score


def fisher(mat):
    se_itr = 0
    adZero = 0
    bcZero = 0
    nonZero = 0
    p_valArr = []
    variance = []
    zAxis = []
    twoXtwoAllSes = [[0, 0], [0, 0]]
    inactive = []
    for se in mat:
        p_valArr2 = []
        variance2 = []
        covArray = []
        enhsLen = len(se["enhs"])
        seRange = se["seId"].split(":")[1].split("-")
        chrNum = se["seId"].split(":")[0].split("chr")[1]
        seSize = int(seRange[1]) - int(seRange[0])
        fibersLen = len(se["enhs"][0]["fibers"])

        twoXtwoSE = [[0, 0], [0, 0]]
        if enhsLen > 1:
            for enh_i in range(enhsLen):
                enh_a = se["enhs"][enh_i]
                thres = 0.10
                inactivePercent1 = len(list(filter(lambda x: x > thres, list(filter(None, enh_a["fibers"]))))) / len(
                    list(filter(None, enh_a["fibers"]))
                )

                column_covariance = []
                for enh_i2 in range(enhsLen):
                    enh_b = se["enhs"][enh_i2]
                    enh_a_range = enh_a["enhId"].split(":")[1].split("-")
                    enh_b_range = enh_b["enhId"].split(":")[1].split("-")
                    enh_dist = int(enh_b_range[0]) - int(enh_a_range[1])
                    inactivePercent2 = len(
                        list(
                            filter(
                                lambda x: x > thres,
                                list(filter(None, enh_b["fibers"])),
                            )
                        )
                    ) / len(list(filter(None, enh_b["fibers"])))

                    if enh_i2 > enh_i and enh_dist < lengthCutOff:
                        inactive.append(inactivePercent1)
                        inactive.append(inactivePercent2)
                        maxA = max(filter(None, enh_a["fibers"]))
                        maxB = max(filter(None, enh_b["fibers"]))
                        # twoXtwo = [[0, 0], [0, 0]]

                        twoXtwo = [[0, 0], [0, 0]]
                        compareA = []
                        compareB = []
                        coverage = 0
                        for fib_i in range(fibersLen):
                            if enh_a["fibers"][fib_i] is not None and enh_b["fibers"][fib_i] is not None:
                                coverage += 1
                                if enh_a["fibers"][fib_i] < thres:
                                    compareA.append(1)
                                else:
                                    compareA.append(0)
                                if enh_b["fibers"][fib_i] < thres:
                                    compareB.append(1)
                                else:
                                    compareB.append(0)
                                if enh_a["fibers"][fib_i] < thres and enh_b["fibers"][fib_i] < thres:
                                    twoXtwo[0][0] += 1
                                    twoXtwoSE[0][0] += 1
                                    twoXtwoAllSes[0][0] += 1
                                if enh_a["fibers"][fib_i] < thres and enh_b["fibers"][fib_i] >= thres:
                                    # bottom left
                                    twoXtwo[1][0] += 1
                                    twoXtwoSE[1][0] += 1
                                    twoXtwoAllSes[1][0] += 1
                                if enh_a["fibers"][fib_i] >= thres and enh_b["fibers"][fib_i] >= thres:
                                    twoXtwo[1][1] += 1
                                    twoXtwoSE[1][1] += 1
                                    twoXtwoAllSes[1][1] += 1
                                if enh_a["fibers"][fib_i] >= thres and enh_b["fibers"][fib_i] < thres:
                                    # top right
                                    twoXtwo[0][1] += 1
                                    twoXtwoSE[0][1] += 1
                                    twoXtwoAllSes[0][1] += 1
                        # jacCov =
                        if coverage > 10:
                            denom = ((twoXtwo[1][0] + twoXtwo[0][1]) / 2) ** 2
                            if denom >= 10:
                                # zAxis.append(jaccard)
                                stat = (twoXtwo[0][0] ** 2) / denom
                                zAxis.append(stat)

    print("      A ON, A OFF")
    print("B ON", twoXtwoAllSes[0])
    print("B OFF", twoXtwoAllSes[1])
    # print("adZeroes: ", adZero, "bcZeroes: ", bcZero, "noZeroes: ", nonZero)
    print("CE Pairs:", len(zAxis))
    print(scipy.stats.mode(zAxis))
    return (lengthCutOff, numpy.median(zAxis), zAxis)


def load_object():
    with open(prefix + "_obj.json") as json_file:
        json_data = json.load(json_file)
        return json_data


def format_matrix(data):
    matrix = []
    enh_ids = []

    for item in data["enhs"]:
        fiber_values = "\t".join("{:<4}".format(str(fiber)) for fiber in item["fibers"])
        print(item["enhId"] + "\t" + fiber_values)


prefixes = [
    "med1_stitched",
    "h3k27ac_Fire",
    "med1_Fire",
    "young_FIRE",
    "youngOurRose_Fire",
]
prefixes = [
    "med1_stitched",
    "med1_Fire",
    "youngOurRose_Fire",
    "youngSE_OSN",
    "youngStitched_OSN",
]
prefixes = [
    # "fire_stitched",
    # "med1_Fire",
    # "youngOurRose_Fire",
    # "youngSE_OSN",
    # "youngStitched_OSN",
]

"""
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
"""
color = {
    "fire_stitched": "#1f77b4",
    "FIRE_peaks_stitched": "#1f77b4",
    "h3k27ac_Fire": "#ff7f0e",
    "med1_Fire": "#d62728",
    "FIRE_peaks_super_med1": "#d62728",
    "young_FIRE": "#9467bd",
    "py_stitch_FIRE_v2": "#9467bd",
    "youngOurRose_Fire": "#2ca02c",
    "young_stitch_FIRE": "#17becf",
    "youngSE_OSN": "#e377c2",
    "youngStitched_OSN": "#8c564b",
    "controls": "#7f7f7f",
}


lengths = [
    10,
    50,
    100,
    200,
    300,
    400,
    500,
    600,
    700,
    800,
    900,
    1000,
    2000,
    3000,
    4000,
    5000,
    6000,
    7000,
    8000,
    9000,
    10000,
    15000,
]


def main_length(prefixes=[], savefig=""):
    mpl_style.set_style()
    global prefix
    global lengthCutOff
    global se_fire_score
    lengthCutOff = 0
    covArrs = []
    plt.figure()
    for prefix in prefixes:
        xAxis = []
        yAxis = []
        yAxisDist = []
        se_fire_score = {}
        mat = load_object()
        # create_super_ranks(mat)
        for lengthCutOff in lengths:
            print(prefix + "_" + str(lengthCutOff))

            xYData = fisher(mat)
            xAxis.append(xYData[0])
            yAxis.append(xYData[1])
            yAxisDist.append(xYData[2])
        if prefix == "controls/FIRE_stitched_control0" or prefix == "osn_control/OSN_stitched_control0":
            prefix = "controls"
            prefix_2 = "Controls"
        fig = plt.figure(figsize=(4, 3))

        ax = fig.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot
        ax.spines[["right", "top"]].set_visible(False)
        ax.spines["bottom"].set_linewidth(1.5)  # Bottom (x-axis) spine
        ax.spines["left"].set_linewidth(1.5)  # Left (y-axis) spine
        ax.tick_params(axis="both", width=2)
        sns.scatterplot(
            x=xAxis,
            y=yAxis,
            # label=prefix,
            color="black",
            alpha=1,
            s=20,
        )
        for x, y, dist in zip(xAxis, yAxis, yAxisDist):
            q25 = np.percentile(dist, 25)
            q75 = np.percentile(dist, 75)
            ax.vlines(x, q25, q75, color="black", linewidth=0.5)
        # print(len(mat))
        # covArrs.append(get_correlation(mat))
        # print(scipy.stats.ks_2samp(covArrs[0],covArrs[1]))
    # plt.legend()
    # plt.legend()
    # l = ax.legend(frameon=True, edgecolor="black", fancybox=False)
    # l.get_frame().set_linewidth(1.0)
    plt.xlabel("Distance between compared enhancers (bp)")
    plt.ylabel("Odds Ratio")
    plt.savefig(savefig)
    plt.show()


main_length(prefixes=["../data/3runsV2"], savefig="../figures/s3b_distance_bin_2575.svg")
# main_length()
# print(scipy.stats.fisher_exact([[2000, 2005], [2000, 1995]]))
