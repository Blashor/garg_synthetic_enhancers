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
import matplotlib

matplotlib.use("Agg")


#
# Calculates co-accessibility statistic for each FIRE peak pair
#


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


def calculate_cutoff(inputVector, drawPlot=True, return_y_thresh=False, savePlot=""):
    inputVector = numpy.sort(inputVector)
    inputVector[inputVector < 0] = 0

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
            f"x={xPt}\ny={y_cutoff:.3f}\nFold over Median={y_cutoff / numpy.median(inputVector):.3f}x\nFold over Mean={y_cutoff / numpy.mean(inputVector):.3f}x"
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


def fisher(mat, dist=False):
    se_itr = 0
    p_valArr = []
    dist_tupple = []
    # has a summary stat as well as identifiers
    summary_ce_pairs = []
    summary_clusters = []

    twoXtwoAllSes = [[0, 0], [0, 0]]
    tupple_arr = []
    thres = 0.10
    linkage_distances = []
    doClusters = False
    drawHeatmaps = False
    writeToFile = True
    for se in mat:
        # print(se["seId"])
        # sox2 "chr3:34716527-34761651"
        # nanog chr6:122686797-122704023
        if se["seId"] != "d":
            other_arr = []
            ce_fisher = []
            ceP_fisher = []
            p_valArr2 = []
            enhsLen = len(se["enhs"])
            seRange = se["seId"].split(":")[1].split("-")
            chrNum = se["seId"].split(":")[0].split("chr")[1]
            seSize = int(seRange[1]) - int(seRange[0])
            fibersLen = len(se["enhs"][0]["fibers"])
            twoXtwoSE = [[0, 0], [0, 0]]
            if enhsLen > 1:
                # print(se["seId"])
                # print(seSize)
                # print(fibersLen)
                fisherAB = [[1 for x in range(enhsLen)] for y in range(enhsLen)]
                fisherAB_2 = [[1 for x in range(enhsLen)] for y in range(enhsLen)]
                fisherAB_3 = [[1 for x in range(enhsLen)] for y in range(enhsLen)]
                se_fisher_mean = []
                # print(se["enhs"])
                for enh_i in range(enhsLen):
                    currentRow = []
                    twoXtwoCE = [[0, 0], [0, 0]]
                    enh_a = se["enhs"][enh_i]
                    enh_a_range = enh_a["enhId"].split(":")[1].split("-")
                    for enh_i2 in range(enhsLen):
                        if enh_i != enh_i2:
                            enh_b = se["enhs"][enh_i2]
                            enh_b_range = enh_b["enhId"].split(":")[1].split("-")
                            if int(enh_a_range[1]) > int(enh_b_range[1]):
                                enh_dist = int(enh_a_range[0]) - int(enh_b_range[1])
                            else:
                                enh_dist = int(enh_b_range[0]) - int(enh_a_range[1])
                            twoXtwo = [[0, 0], [0, 0]]
                            enh_a_state = []
                            enh_b_state = []
                            for fib_i in range(fibersLen):
                                if enh_a["fibers"][fib_i] is not None and enh_b["fibers"][fib_i] is not None:
                                    if enh_a["fibers"][fib_i] < thres:
                                        enh_a_state.append(1)
                                    else:
                                        enh_a_state.append(0)
                                    if enh_b["fibers"][fib_i] < thres:
                                        enh_b_state.append(1)
                                    else:
                                        enh_b_state.append(0)
                                    if enh_a["fibers"][fib_i] < thres and enh_b["fibers"][fib_i] < thres:
                                        twoXtwo[0][0] += 1
                                        twoXtwoSE[0][0] += 1
                                        twoXtwoAllSes[0][0] += 1
                                        twoXtwoCE[0][0] += 1
                                    if enh_a["fibers"][fib_i] < thres and enh_b["fibers"][fib_i] >= thres:
                                        # bottom left
                                        twoXtwo[1][0] += 1
                                        twoXtwoSE[1][0] += 1
                                        twoXtwoAllSes[1][0] += 1
                                        twoXtwoCE[1][0] += 1
                                    if enh_a["fibers"][fib_i] >= thres and enh_b["fibers"][fib_i] >= thres:
                                        twoXtwo[1][1] += 1
                                        twoXtwoSE[1][1] += 1
                                        twoXtwoAllSes[1][1] += 1
                                        twoXtwoCE[1][1] += 1
                                    if enh_a["fibers"][fib_i] >= thres and enh_b["fibers"][fib_i] < thres:
                                        # top right
                                        twoXtwo[0][1] += 1
                                        twoXtwoSE[0][1] += 1
                                        twoXtwoAllSes[0][1] += 1
                                        twoXtwoCE[0][1] += 1
                            a = twoXtwo[0][0]
                            b = twoXtwo[0][1]
                            c = twoXtwo[1][0]
                            d = twoXtwo[1][1]
                            n = a + b + c + d
                            # print(a, b, c, d)
                            if a + b + c + d >= 10 and b + c > 0 and a > 0:
                                # binomResult = scipy.stats.binomtest(a, n, p=1 / 3, alternative="greater")
                                pvalue = 1  # binomResult.pvalue

                                # print()
                                # print(enh_i, enh_i2, a, b, c, binomResult.statistic)
                                # unbounded comparison
                                jaccard_sym = (2 * a) / (2 * a + b + c)

                                """concord_discord = (
                                    4 * numpy.cov(enh_a_state, enh_b_state)[0][1]
                                )
                                """
                                concord_discord = (a) ** 2 / ((b + c) / 2) ** 2
                                if isinstance(dist, tuple):
                                    if enh_dist > 50:
                                        expected = powlaw(enh_dist, *dist)
                                        concord_discord = concord_discord / expected
                                    else:
                                        concord_discord = 0
                                # concord_discord = pearsonr(enh_a_state, enh_b_state)[0]
                                other_score = (str(a), str(b), str(c), str(d))
                                #
                                #
                                #
                                p_score = -10 * math.log10(pvalue)
                                # concord_discord = p_score
                                # concord_discord = (2 * a) / (b + c)
                                # concord_discord = a / (a + b + c)
                                # print(enh_dist, concord_discord)
                                if dist == True:
                                    dist_tupple.append((enh_dist, concord_discord))
                                summary_ce_pairs.append(
                                    (
                                        concord_discord,
                                        se["seId"],
                                        enh_a["enhId"],
                                        enh_b["enhId"],
                                    )
                                )
                                # jaccard = binomResult.statistic
                                # symmetrical version purely for euclidean clustering

                                if pvalue < 1.1:
                                    fisherAB[enh_i][enh_i2] = jaccard_sym
                                    if drawHeatmaps:
                                        print(
                                            enh_i,
                                            enh_i2,
                                            concord_discord,
                                            a,
                                            b,
                                            c,
                                            d,
                                        )

                                    fisherAB_2[enh_i][enh_i2] = concord_discord
                                    fisherAB_3[enh_i][enh_i2] = other_score
                                    se_fisher_mean.append(jaccard_sym)
                                else:
                                    fisherAB[enh_i][enh_i2] = 0.5
                            else:
                                fisherAB[enh_i][enh_i2] = 0.5
                        else:
                            fisherAB[enh_i][enh_i2] = 0
                    twoXtwoCE[0][1] = (twoXtwoCE[0][1] + twoXtwoCE[1][0]) / 2
                    twoXtwoCE[1][0] = twoXtwoCE[0][1]
                    twoXtwoCE[1][1] = twoXtwoCE[0][0]
                    statisticCE, p_valCE = scipy.stats.fisher_exact(twoXtwoCE)
                    ce_fisher.append(round(statisticCE, 2))
                    ceP_fisher.append(round(p_valCE, 2))
                    if p_valCE < 0.25:
                        other_arr.append(statisticCE)
                    else:
                        pass
                newFisher = []

                for row in fisherAB:
                    localMax = numpy.max(row)
                    row = numpy.where(numpy.isclose(row, 0), localMax, row)
                    newFisher.append(row)
                fisherAB = newFisher
                se_mean = str(round(numpy.mean(se_fisher_mean), 2))
                mask = numpy.triu(fisherAB)
                seperate_linkage = fastcluster.linkage_vector(fisherAB, metric="euclidean", method="ward")
                if doClusters == True and dist == False:
                    # print(seperate_linkage)
                    cut_off = get_dendogram_cutoff(seperate_linkage)

                    # Clusters node by distance
                    clusters = []
                    nodes_to_explore = []
                    rootnode = scipy.cluster.hierarchy.to_tree(seperate_linkage)
                    nodes_to_explore.append(rootnode)
                    while len(nodes_to_explore) > 0:
                        node = nodes_to_explore[0]
                        right_child = node.right
                        left_child = node.left

                        if node.dist >= cut_off:
                            # print(node.id)
                            if right_child is not None:
                                nodes_to_explore.append(right_child)
                            if left_child is not None:
                                nodes_to_explore.append(left_child)
                        else:
                            # print("cluster", str(node.id))
                            clusters.append((node.pre_order(lambda x: x.id)))
                        nodes_to_explore = nodes_to_explore[1:]
                    # print(clusters)
                else:
                    clusters = [list(range(len(fisherAB)))]
                for cluster in clusters:
                    ce_string = ""
                    cluster_scores = []
                    other_scores = []
                    for enh_i in cluster:
                        ce_string += se["enhs"][enh_i]["enhId"] + ", "
                        for enh_i2 in cluster:
                            if enh_i < enh_i2:
                                # print(enh_i, enh_i2, fisherAB[enh_i][enh_i2])
                                if fisherAB[enh_i][enh_i2] != 0.5:
                                    cluster_scores.append(fisherAB_2[enh_i][enh_i2])
                                    other_scores.append(fisherAB_3[enh_i][enh_i2])
                                    """
                                    print(
                                        fisherAB_2[enh_i][enh_i2],
                                        fisherAB[enh_i][enh_i2],
                                        enh_i,
                                        enh_i2,
                                    )
                                    """
                    if len(cluster_scores) > 0:
                        # c_mean = numpy.mean(cluster_scores)
                        c_max = numpy.max(cluster_scores)

                        max_itr = numpy.array(cluster_scores).argmax()

                        summary_clusters.append((c_max, se["seId"], ce_string, other_scores[max_itr]))
                        # print(c_mean)

                if drawHeatmaps == True:
                    col_map = sns.color_palette()
                    color_dict = {}
                    # row_colors = pd.DataFrame()[0].map(color_dict)
                    for c_itr in range(len(clusters)):
                        cluster = clusters[c_itr]
                        color = col_map[c_itr]
                        for leaf in cluster:
                            color_dict[leaf] = color
                    list_num = list(range(len(fisherAB)))
                    df = pd.DataFrame(fisherAB_2, index=list_num)
                    row_colors = pd.DataFrame(list_num)[0].map(color_dict)
                    #
                    cluster_grid = sns.clustermap(
                        df,
                        annot=True,
                        cmap="seismic",
                        vmin=0,
                        vmax=36,
                        center=1,
                        row_colors=row_colors,
                        col_colors=row_colors,
                        row_linkage=seperate_linkage,
                        col_linkage=seperate_linkage,
                        linewidths=1,
                        # tree_kws={"colors": [col_map[s - 1] for s in cut_tree]},
                        linecolor="white",
                        # method="ward",
                        figsize=(12, 6),
                        # ax=axs[0]
                    )
                    cluster_grid.fig.subplots_adjust(right=0.5)
                    ax = cluster_grid.fig.add_axes([0.6, 0.05, 0.42, 0.8])
                    sns.heatmap(
                        fisherAB_2,
                        square=True,
                        cmap="seismic",
                        center=1,
                        vmax=36,
                        annot=True,
                        linewidths=1,
                        linecolor="white",
                        mask=mask,
                        ax=ax,
                        cbar=True,
                    )
                    cluster_grid.fig.suptitle(se["seId"])
                    plt.yticks(rotation=0)
                    plt.show()

                se_itr += 1

    if dist == True:
        # random.Random(0).shuffle(dist_tupple)
        dist_tupple.sort(key=lambda tup: tup[0])
        return dist_tupple
    #
    # SAVING CEs
    #

    if writeToFile == True:
        summary_ce_pairs.sort(key=lambda tup: tup[0], reverse=True)
        ce_tupples = list(map(list, zip(*summary_ce_pairs)))
        y_cutoff = calculate_cutoff(
            ce_tupples[0],
            drawPlot=True,
            return_y_thresh=True,
            savePlot="data/ce_pairs_ranked.svg",
        )
        ce_lines = []
        rankNum = 0
        for ce_pair in summary_ce_pairs:
            isSuper = ""
            rankNum += 1
            if ce_pair[0] >= y_cutoff:
                isSuper = "Super"
            ce_lines.append(
                ce_pair[2]
                + "\t"
                + ce_pair[3]
                + "\t"
                + ce_pair[1]
                + "\t"
                + str(ce_pair[0])
                + "\t"
                + str(rankNum)
                + "\t"
                + isSuper
                + "\n"
            )

        with open("data/ce_rank.txt", "w") as file:
            file.writelines(ce_lines)

        #
        # SAVING Clusters
        #
        summary_clusters.sort(key=lambda tup: tup[0], reverse=True)
        cluster_tupples = list(map(list, zip(*summary_clusters)))
        y_cutoff = calculate_cutoff(
            cluster_tupples[0],
            drawPlot=True,
            return_y_thresh=True,
            savePlot="data/clusters_ranked.svg",
        )
        cluster_lines = []
        rankNum = 0
        for cluster in summary_clusters:
            isSuper = ""
            rankNum += 1
            constituent_num = len(cluster[2].split(",")) - 1
            if cluster[0] >= y_cutoff:
                isSuper = "Super"
            cluster_lines.append(
                cluster[2]
                + "\t"
                + cluster[1]
                + "\t"
                + str(constituent_num)
                + "\t"
                + str(cluster[0])
                + "\t"
                + str(rankNum)
                + "\t"
                + isSuper
                + "\t"
                + "\t".join(cluster[3])
                + "\n"
            )

        with open("data/cluster_rank.txt", "w") as file:
            file.writelines(cluster_lines)

    print(prefix)
    print("      A ON, A OFF")
    print("B ON", twoXtwoAllSes[0])
    print("B OFF", twoXtwoAllSes[1])
    fishResAllSes = scipy.stats.fisher_exact(twoXtwoAllSes)
    print(fishResAllSes)
    print("SEs: " + str(len(p_valArr)))
    sns.ecdfplot(data=linkage_distances)
    plt.show()
    # tupple_arr.sort(key=lambda tup: tup[2], reverse=False)
    # print(tupple_arr)
    """
    set1 = set(super_arr)
    set2 = set(super_arr_2)
    venn2([set1, set2], set_labels=("CE", "SE"))

    plt.show()
    """


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


def load_object():
    with open(prefix + "_obj.json") as json_file:
        json_data = json.load(json_file)
        return json_data


prefixes = ["/Users/blake/Downloads/K562"]
# prefixes = ["FIRE_peaks_stitch_ROSElike"]


def main():
    global prefix
    # prefixes = sys.argv[1:]
    for prefix in prefixes:
        print(prefix)
        mat = load_object()
        # fisher(mat)
        dist_tupple = fisher(mat, dist=True)
        print(dist_tupple)
        dist_correction = mcfs_dist(dist_tupple)
        fisher(mat, dist=dist_correction)


# print(sns.cm.rocket)
main()
