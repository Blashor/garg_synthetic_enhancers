import math
import sys
import json
import scipy
import numpy
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2
import fastcluster
from sklearn.cluster import KMeans
import pandas as pd
from scipy.stats.stats import pearsonr
import random
from scipy.optimize import curve_fit


def get_dendogram_cutoff(data):
    if len(data) == 1:
        return 1
    node_lookup = []
    node_dist = []
    nodes_to_explore = []
    rootnode = scipy.cluster.hierarchy.to_tree(data)
    nodes_to_explore.append(rootnode)
    while len(nodes_to_explore) > 0:
        node = nodes_to_explore[0]
        right_child = node.right
        left_child = node.left
        dist_jump = []
        child_num = 0
        if right_child is not None:
            dist_jump.append(node.dist - right_child.dist)

            nodes_to_explore.append(right_child)
        if left_child is not None:
            dist_jump.append(node.dist - left_child.dist)

            nodes_to_explore.append(left_child)
        if len(dist_jump) > 0:
            node_dist.append(numpy.mean(dist_jump))
            node_lookup.append(node.dist)
        nodes_to_explore = nodes_to_explore[1:]
    max_jump_node = node_lookup[numpy.argmax(numpy.array(node_dist))]
    pearson_r = scipy.stats.pearsonr(range(len(data[:, 2])), data[:, 2]).statistic
    # print(max_jump_node)
    # print(pearson_r)
    if pearson_r > 0.99:
        return max(data[:, 2]) + 1
    return max_jump_node


def numPts_below_line(myVector, slope, x):
    yPt = myVector[x]
    b = yPt - (slope * x)
    xPts = numpy.arange(1, len(myVector) + 1)
    return numpy.sum(myVector <= (xPts * slope + b))


def calculate_cutoff(inputVector, label, drawPlot=True, return_y_thresh=False, savePlot=""):
    inputVector = numpy.sort(inputVector)
    inputVector[inputVector < 0] = 0
    n_len = len(inputVector)
    slope = (numpy.max(inputVector) - numpy.min(inputVector)) / len(inputVector)

    result = numpy.vectorize(lambda x: numPts_below_line(inputVector, slope, x))(numpy.arange(len(inputVector)))
    xPt = numpy.floor(numpy.argmin(result)) + 1
    y_cutoff = inputVector[int(xPt - 1)]

    if drawPlot:
        plt.figure(figsize=(10, 6))
        plt.plot(numpy.arange(1, len(inputVector) + 1), inputVector, "b-")
        b = y_cutoff - (slope * xPt)
        # plt.axvline(x=xPt, color="gray", linestyle="dashed")
        # plt.axhline(y=y_cutoff, color="gray", linestyle="dashed")
        plt.scatter(xPt, y_cutoff, marker="o", color="red")
        """
        plt.plot(
            numpy.arange(1, len(inputVector) + 1),
            slope * numpy.arange(1, len(inputVector) + 1) + b,
            "r-",
        )
        """
        plt.title(
            f"{label}\nn={n_len}, x={xPt}, y={y_cutoff:.3f}\nFold over Median={y_cutoff / numpy.median(inputVector):.3f}x\nFold over Mean={y_cutoff / numpy.mean(inputVector):.3f}x"
        )
        """
        plt.axvline(
            x=numpy.sum(inputVector == 0), color="pink", linestyle="dotted", linewidth=2
        )
        """
        plt.xlabel("Enhancers Ranked By Modified Odds Ratio")
        plt.ylabel("Modified Odds Ratio")
        if savePlot != "":
            plt.savefig(savePlot)
    if return_y_thresh == True:
        return y_cutoff
    return xPt


def fisher(mat, label, dist=False):
    se_itr = 0
    thres = 0.10
    ce_fisher_scores = []
    ce_output = []
    for se in mat:
        if se["seId"] != "_chr3:34716527-34761651":
            enhsLen = len(se["enhs"])
            fibersLen = len(se["enhs"][0]["fibers"])
            if enhsLen > 1:
                for enh_i in range(enhsLen):
                    twoXtwoCE = [[0, 0], [0, 0]]
                    enh_a = se["enhs"][enh_i]
                    enh_a_id = enh_a["enhId"]
                    a_start, a_end = map(int, enh_a_id.split(":")[1].split("-"))
                    for enh_i2 in range(enhsLen):
                        if enh_i < enh_i2:
                            enh_b = se["enhs"][enh_i2]
                            enh_b_id = enh_b["enhId"]
                            b_start, b_end = map(int, enh_b_id.split(":")[1].split("-"))
                            twoXtwo = [[0, 0], [0, 0]]
                            for fib_i in range(fibersLen):
                                if enh_a["fibers"][fib_i] is not None and enh_b["fibers"][fib_i] is not None:
                                    if enh_a["fibers"][fib_i] < thres and enh_b["fibers"][fib_i] < thres:
                                        twoXtwo[0][0] += 1
                                    if enh_a["fibers"][fib_i] < thres and enh_b["fibers"][fib_i] >= thres:
                                        # bottom left
                                        twoXtwo[1][0] += 1
                                    if enh_a["fibers"][fib_i] >= thres and enh_b["fibers"][fib_i] >= thres:
                                        twoXtwo[1][1] += 1
                                    if enh_a["fibers"][fib_i] >= thres and enh_b["fibers"][fib_i] < thres:
                                        # top right
                                        twoXtwo[0][1] += 1
                            a = twoXtwo[0][0]
                            b = twoXtwo[0][1]
                            c = twoXtwo[1][0]
                            d = twoXtwo[1][1]
                            n = a + b + c + d
                            if n > 10:
                                denom = ((b + c) * 0.5) ** 2
                                if denom != 0:
                                    bs_stat = (a**2) / denom
                                    ce_fisher_scores.append(bs_stat)
                                    abcd = "\t".join(map(str, ([a, b, c, d])))
                                    out_line = f"{enh_a_id}\t{enh_b_id}\t{label}\t{abcd}\t{bs_stat}"
                                    ce_output.append(out_line + "\n")
    print(label, len(ce_fisher_scores))
    with open(f"data/{label}_ce_score.txt", "w") as file:
        file.writelines(ce_output)
    # calculate_cutoff(ce_fisher_scores, label)
    sns.ecdfplot(data=ce_fisher_scores, label=label)


def powlaw(x, a, b, c):
    return a / numpy.power(x, b) + c


def mcfs_dist(dist_tupple):
    window_size = 50
    windows_span = 0
    median_arr = []
    dist_arr = []
    x = []
    y = []
    for pair in dist_tupple:
        dist_arr.append(pair[0])
        median_arr.append(pair[1])
        if len(median_arr) == 500:
            print(numpy.median(dist_arr), numpy.percentile(median_arr, [25, 50, 75]))
            # print(median_arr)
            x.append(numpy.median(dist_arr))
            y.append(numpy.median(median_arr))
            median_arr = []
            dist_arr = []
    new_x = numpy.linspace(x[0], x[-1], int(0.5 * len(x)))
    new_y = numpy.interp(new_x, x, y)
    sns.scatterplot(x=x, y=y)
    # sns.scatterplot(x=new_x, y=new_y, color=(0, 0, 0, 0.2))
    sol = curve_fit(powlaw, new_x, new_y, maxfev=10000000)
    print(sol)
    a, b, c = sol[0]
    y = powlaw(x, a, b, c)
    sns.lineplot(x=x, y=y)
    plt.show()
    return (a, b, c)


def cluster_mat(filename, mat):
    cluster_mats = {}
    cluster_labels = {}
    with open(filename) as file:
        for line in file:
            if "deepTools_group" not in line:
                cell = line.strip().split("\t")
                ce_name = f"{cell[0]}:{cell[1]}-{cell[2]}"
                cluster_labels[ce_name] = cell[-1]
                if cell[-1] not in cluster_mats:
                    cluster_mats[cell[-1]] = []
    for stitch in mat:
        temp_cluster_holder = {}
        for ce in stitch["enhs"]:
            ce_cluster = cluster_labels[ce["enhId"]]
            if ce_cluster not in temp_cluster_holder:
                temp_cluster_holder[ce_cluster] = {
                    "chr": ce["enhId"].split(":")[0],
                    "min": 10**12,
                    "max": 0,
                    "enhs": [],
                }
                start, end = map(int, ce["enhId"].split(":")[1].split("-"))
                if start < temp_cluster_holder[ce_cluster]["min"]:
                    temp_cluster_holder[ce_cluster]["min"] = start
                if end > temp_cluster_holder[ce_cluster]["max"]:
                    temp_cluster_holder[ce_cluster]["max"] = end
            temp_cluster_holder[ce_cluster]["enhs"].append(ce)
        for cluster in temp_cluster_holder:
            c_handle = temp_cluster_holder[cluster]
            chrom = c_handle["chr"]
            start = c_handle["min"]
            end = c_handle["max"]
            cluster_se_id = f"{chrom}:{start}-{end}"
            new_stitch_obj = {"seId": cluster_se_id, "enhs": c_handle["enhs"]}
            cluster_mats[cluster].append(new_stitch_obj)

    return cluster_mats


def load_object():
    with open(prefix + "_obj.json") as json_file:
        json_data = json.load(json_file)
        return json_data


prefixes = ["../FIRE/FIRE_peaks_stitched_v2"]
# prefixes = ["FIRE_peaks_stitch_ROSElike"]


def main():
    global prefix
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot
    ax.spines[["right", "top"]].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.5)  # Bottom (x-axis) spine
    ax.spines["left"].set_linewidth(1.5)  # Left (y-axis) spine
    ax.tick_params(axis="both", width=2)
    for prefix in prefixes:
        print(prefix)
        mat = load_object()
        clusters = cluster_mat("../metagene/fire_ce.svg3.bed", mat)
        for cl in clusters:
            fisher(clusters[cl], cl)
        fisher(mat, "Unclustered")
        plt.legend()
        plt.xscale("log")
        plt.axvline(x=1, color="k", linestyle="--")

        plt.show()


# print(sns.cm.rocket)
main()
