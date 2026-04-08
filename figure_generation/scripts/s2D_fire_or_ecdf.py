import math
import sys
import json
import scipy
import numpy
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np


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
    coverageArr = []
    coverageArrAfterFilt = []
    p_valArr = []
    data_stats_arr = [[], [], [], []]
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
        # twoXtwoSE = [[1, 1], [1, 1]]
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

                    if enh_i2 > enh_i:
                        inactive.append(inactivePercent1)
                        inactive.append(inactivePercent2)
                        maxA = max(filter(None, enh_a["fibers"]))
                        maxB = max(filter(None, enh_b["fibers"]))
                        # twoXtwo = [[0, 0], [0, 0]]
                        # print(se)
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
                                    twoXtwoAllSes[1][0] += 1
                                    twoXtwoSE[1][0] += 1
                                if enh_a["fibers"][fib_i] >= thres and enh_b["fibers"][fib_i] >= thres:
                                    twoXtwo[1][1] += 1
                                    twoXtwoAllSes[1][1] += 1
                                    twoXtwoSE[1][1] += 1
                                if enh_a["fibers"][fib_i] >= thres and enh_b["fibers"][fib_i] < thres:
                                    # top right
                                    twoXtwo[0][1] += 1
                                    twoXtwoAllSes[0][1] += 1
                                    twoXtwoSE[0][1] += 1
                        coverageArr.append(coverage)

                        # enhancer level adding to graph
                        '''
                        fishRes = scipy.stats.fisher_exact(twoXtwo)
                        p_val = fishRes[1]
                        statistic = fishRes[0]
                        if str(statistic) != "inf" and str(statistic) != "nan":
                            if twoXtwoSE[0][0] < twoXtwoSE[1][1]:
                                pass
                                # statistic *= -1
                            """var = numpy.var(
                                numpy.ndarray.flatten(numpy.array(twoXtwoSE))
                            )"""
                            # cov = numpy.cov(m=compareA, y=compareB)[1][1]
                            if twoXtwo[1][1] != 0 and twoXtwo[0][0] != 0:
                                onOffRatio = twoXtwo[0][0] / twoXtwo[1][1]
                                # variance.append(math.log(var + 1))
                                if coverage > 4:
                                    # zAxis.append(enh_dist)
                                    variance.append(onOffRatio)
                                    p_valArr.append(statistic)
                        '''
            se_itr += 1
            # print(twoXtwoSE)

            if twoXtwoSE[1][1] == 0 or twoXtwoSE[0][0] == 0:
                adZero += 1
            elif twoXtwoSE[0][1] == 0 or twoXtwoSE[1][0] == 0:
                bcZero += 1
            else:
                nonZero += 1
            # twoXtwoSE[1][1] = twoXtwoSE[0][0]
            # print(twoXtwoSE)
            # twoXtwoSE[1][0] = int((twoXtwoSE[1][0] + twoXtwoSE[0][1]) / 2)
            # twoXtwoSE[0][1] = twoXtwoSE[1][0]
            # print(twoXtwoSE)
            statistic, p_val = scipy.stats.fisher_exact(twoXtwoSE)
            if str(statistic) != "inf" and str(statistic) != "nan":
                if twoXtwoSE[1][1] != 0 and twoXtwoSE[0][0] != 0:
                    a = twoXtwoSE[0][0]
                    b = twoXtwoSE[0][1]
                    c = twoXtwoSE[1][0]
                    bc = twoXtwoSE[1][0] * twoXtwoSE[0][1]
                    d = twoXtwoSE[1][1]
                    w = 0.25
                    q = 1 / (1 - w**2)
                    avg = (twoXtwoSE[0][1] + twoXtwoSE[1][0]) / 2
                    a2OverBC = (twoXtwoSE[0][0] ** 2) / (avg) ** 2
                    aOverBplusC = twoXtwoSE[0][0] / (twoXtwoSE[0][1] + twoXtwoSE[1][0])
                    jaccard = twoXtwoSE[0][0] / (twoXtwoSE[0][1] + twoXtwoSE[1][0] + twoXtwoSE[0][0])

                    # print(addRatio, twoXtwoSE)
                    if p_val < 1.5 and a + b + c > 10:
                        data_stats_arr[0].append(0)
                        data_stats_arr[1].append(a2OverBC)
                        data_stats_arr[2].append(0)
                        data_stats_arr[3].append(0)
                    # p_valArr.append(a2OverBC)
                # p_valArr.append(onOffRatio)

            #    se_itr += 1
            # sys.exit()
    binsNum = 250
    # print(len(p_valArr))
    # print(p_valArr)
    # twoXtwoAllSes[1][1] = int((twoXtwoAllSes[0][1] + twoXtwoAllSes[1][0]) / 2)

    print("      A ON, A OFF")
    print("B ON", twoXtwoAllSes[0])
    print("B OFF", twoXtwoAllSes[1])
    print("adZeroes: ", adZero, "bcZeroes: ", bcZero, "noZeroes: ", nonZero)
    print(
        "Raw Coverage:",
        round(numpy.mean(coverageArr), 2),
        "| Filtered Coverage:",
        round(numpy.mean(coverageArrAfterFilt), 2),
    )
    fishResAllSes = scipy.stats.fisher_exact(twoXtwoAllSes)
    # odds = scipy.stats.contingency.odds_ratio(twoXtwoAllSes)
    # print(odds.confidence_interval(confidence_level=0.95))

    print(fishResAllSes)
    print("SEs: " + str(len(p_valArr)))

    if prefix == "youngSE_OSN" or prefix == "med1_Fire":
        prefix_2 = "Super"
    elif prefix == "youngStitched_OSN" or prefix == "fire_stitched":
        prefix_2 = "Stitched"
    # for itr in range(4):
    itr = 1
    # plt.figure(itr + 1)

    if doingControls == True:
        thePlot = sns.ecdfplot(data=jitter_array(data_stats_arr[itr]), alpha=1, color="#7f7f7f", label="Controls")
    else:
        """
        supers = []
        y_cut = calculate_cutoff(data_stats_arr[itr])
        for val in data_stats_arr[itr]:
            if val >= y_cut:
                supers.append(val)
        thePlot = sns.ecdfplot(
            data=supers, alpha=1, label="Super", color=color["med1_Fire"]
        )
        """
        thePlot = sns.ecdfplot(data=jitter_array(data_stats_arr[itr]), alpha=1, label=prefix_2, color=color[prefix])

        plt.xscale("log")
        plt.legend()
        if itr == 0:
            plt.xlabel("Odds Ratio")
        elif itr == 1:
            plt.xlabel("Odds Ratio")
        thePlot.set_ylim(0, 1.1)
        thePlot.set_xlim(10 ** (-3), 10**3)

    """
    if prefix == "med1_stitched":
        sns.scatterplot(
            x=variance,
            y=p_valArr,
            label=prefix + "_" + str(lengthCutOff),
            alpha=0.2,
            s=10,
        )
    else:
        sns.scatterplot(
            x=variance,
            y=p_valArr,
            label=prefix + "_" + str(lengthCutOff),
            alpha=0.8,
            s=10,
        )
        """


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


# prefix = "H3k27ac_Fire"
# prefix = "youngSE_Fire"


# prefix = "med1_Fire"

# prefixes = ["stitchedControl", "med1_Fire", "H3k27ac_Fire", "youngSE_Fire"]

prefixes = ["med1_stitched", "med1_Fire", "youngOurRose_Fire"]

prefixes = [
    "med1_stitched",
    "h3k27ac_Fire",
    "youngOurRose_Fire",
    "med1_Fire",
    "young_FIRE",
    "youngSE_OSN",
    "youngStitched_OSN",
    "exons_stitched",
    "only_intergenic_stitched",
    "inverse_stitched",
]


prefixes = [
    "med1_stitched",
    "med1_Fire",
    "youngOurRose_Fire",
    "young_stitch_FIRE",
    "youngSE_OSN",
    "youngStitched_OSN",
]


prefixes = [
    "young_FIRE",
    "young_stitch_FIRE",
    "youngSE_OSN",
    "youngStitched_OSN",
]
prefixes = ["fire_stitched", "young_stitch_FIRE", "youngStitched_OSN"]

color = {
    "fire_stitched": "#1f77b4",
    "h3k27ac_Fire": "#ff7f0e",
    "med1_Fire": "#d62728",
    "young_FIRE": "#9467bd",
    "py_stitch_FIRE": "#9467bd",
    "youngOurRose_Fire": "#2ca02c",
    "young_stitch_FIRE": "#17becf",
    "youngSE_OSN": "#e377c2",
    "youngStitched_OSN": "#8c564b",
    "control": "#7f7f7f",
}


prefixes = [
    "youngSE_OSN",
    "youngStitched_OSN",
]

prefixes = ["fire_stitched", "med1_Fire", "py_stitch_FIRE"]


prefixes = [
    "youngSE_OSN",
    "youngStitched_OSN",
]
prefixes = ["fire_stitched", "med1_Fire"]
prefixes = ["fire_stitched"]
prefixes = [
    "youngStitched_OSN",
    "youngSE_OSN",
]
prefixes = ["fire_stitched", "med1_Fire"]

prefixes = ["fire_stitched"]


def jitter_array(arr, jitter_strength=0.03):
    arr = np.array(arr)
    log_arr = np.log10(arr + 1e-8)  # avoid log(0)
    jitter = np.random.normal(loc=0, scale=jitter_strength, size=arr.shape)
    log_arr_jittered = log_arr + jitter
    jittered_arr = 10**log_arr_jittered
    return np.clip(jittered_arr, a_min=1e-3, a_max=1e3)


def main():
    global prefix
    global lengthCutOff
    global se_fire_score
    global doingControls
    lengthCutOff = 0
    covArrs = []
    doingControls = True
    fig = plt.figure(figsize=(4, 3))

    ax = fig.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot
    ax.spines[["right", "top"]].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.5)  # Bottom (x-axis) spine
    ax.spines["left"].set_linewidth(1.5)  # Left (y-axis) spine
    ax.tick_params(axis="both", width=2)
    for itr in range(1):
        print("control: ", itr)
        # prefix = "osn_control/OSN_stitched_control" + str(itr)
        prefix = "controls/FIRE_stitched_control" + str(itr)
        mat = load_object()
        fisher(mat)
    doingControls = False
    itr = 0
    for prefix in prefixes:
        print(prefix, "\n")
        se_fire_score = {}
        # format_matrix(mat[303])
        # Object to correlate, doSpearmanRanking?
        # mat = prep_object()
        mat = load_object()
        fisher(mat)

    # plt.xlabel("log(Odds Ratio)")

    plt.show()
    plt.savefig("../figures/control_odds_ratio.svg")


main()
