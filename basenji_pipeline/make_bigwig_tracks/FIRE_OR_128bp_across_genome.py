# workflow
# 1. binsToBed
# 2. bedtools intersect
# 3. create a matrix from the fiber level of intersection
# 4. dump to pre

import scipy
from scipy.optimize import curve_fit
import matplotlib as mpl
import numpy
import sys

mpl.use("Agg")
from matplotlib import pyplot as plt
import seaborn as sns

from multiprocessing.pool import ThreadPool
import threading
import json

file_lock = threading.Lock()


def binsToBed(bin_size=128):
    bed_lines = []
    #
    with open("/home/bmt26/scripts/CHIPseq/exclusionGFF/mm10Sort.genome") as file:
        for line in file:
            cell = line.split("\t")
            if len(cell[0].split("_")) == 1:
                print(cell)
                for window in range(bin_size, int(cell[1]), bin_size):
                    bed_lines.append("\t".join((cell[0], str(window - bin_size), str(window))) + "\n")
                bed_lines.append("\t".join((cell[0], str(window), cell[1])))
    #
    with open("/home/bmt26/garg/fiberSeq/topology/mm10_" + str(bin_size) + "bp_bins.bed", "w") as file:
        file.writelines(bed_lines)
        print(len(bed_lines))


def to_json(obj, fileName):
    with open(fileName, "w") as file:
        file.write(json.dumps(obj))


def FIRE_intersect_to_matrix(
    topology_path="/home/bmt26/garg/fiberSeq/topology",
    input_bed="bin128_fire_intersect.bed",
    output_short="bin128_whole_genome_fire_mm10.short",
):
    intersect_obj = {}
    old_chrom = ""
    output = f"{topology_path}/"
    with open(output, "w") as out_file:
        pass
    with ThreadPool(2) as pool:
        with open(f"{topology_path}/{input_bed}") as file:
            for line in file:
                cell = line.split("\t")
                # skip to the part
                if len(cell) == 15:
                    # this actually hppens after the json contents are made, so code below executes before "if" ever time
                    if old_chrom != cell[0]:
                        # print(intersect_obj)
                        json_file = f"{topology_path}/" + input_bed.replace(".bed", "") + old_chrom + ".json"
                        to_json(intersect_obj, json_file)
                        pool.apply_async(abcd_matrix_to_file, args=(json_file, output))
                        old_chrom = cell[0]
                        intersect_obj = {}
                    fiber_id = cell[6]
                    centered_coord = int((int(cell[1]) + int(cell[2])) / 2)
                    window = cell[0] + ":" + str(centered_coord)
                    FIRE_score = float(cell[12])
                    if fiber_id not in intersect_obj:
                        intersect_obj[fiber_id] = {}
                    if window not in intersect_obj[fiber_id]:
                        if FIRE_score < 0.1:
                            intersect_obj[fiber_id][window] = int(cell[4]) + int(cell[5])
                        else:
                            intersect_obj[fiber_id][window] = 0
                    else:
                        if FIRE_score < 0.1:
                            peak_id = int(cell[4]) + int(cell[5])
                            if intersect_obj[fiber_id][window] == 0:
                                intersect_obj[fiber_id][window] = peak_id
                            elif isinstance(intersect_obj[fiber_id][window], int):
                                intersect_obj[fiber_id][window] = [
                                    intersect_obj[fiber_id][window],
                                    peak_id,
                                ]
                            else:
                                intersect_obj[fiber_id][window].append(peak_id)

        json_file = f"{topology_path}/" + input_bed.replace(".bed", "") + "_last_chrom.json"
        to_json(intersect_obj, json_file)
        pool.apply_async(abcd_matrix_to_file, args=(json_file, output))
        pool.close()
        pool.join()
        sys.stderr.flush()
        sys.stdout.flush()


def load_object(file_name):
    with open(file_name) as json_file:
        return json.load(json_file)


def abcd_matrix_to_file(json_name, output):
    print(json_name)
    sys.stdout.flush()
    intersect_obj = load_object(json_name)
    sparse_matrix = {}
    out_lines = []
    for fiber_name in intersect_obj:
        intersect_obj[fiber_name] = dict(
            sorted(
                intersect_obj[fiber_name].items(),
                key=lambda item: (item[1] == 0, item[1] == 0),
            )
        )
        fiber = intersect_obj[fiber_name]
        first_time_in_fiber = True
        for window1_key in fiber:
            window1 = fiber[window1_key]
            if first_time_in_fiber == True:
                first_time_in_fiber = False
                # indicates that this fiber has no fire elements
                if window1 == 0:
                    break
            for window2_key in fiber:
                window2 = fiber[window2_key]
                if window1_key != window2_key:
                    if window1 != 0:
                        bin_compare = window1_key + ":" + window2_key
                        if bin_compare not in sparse_matrix:
                            sparse_matrix[bin_compare] = [0, 0, 0, True]
                        elif sparse_matrix[bin_compare][3] == False:
                            continue
                        if window2 != 0:
                            # a
                            if window1 != window2:
                                sparse_matrix[bin_compare][0] += 1
                            else:
                                # blacklisted
                                sparse_matrix[bin_compare][3] = False
                        else:
                            # c
                            sparse_matrix[bin_compare][2] += 1
                    elif window2 != 0:
                        # b
                        bin_compare = window1_key + ":" + window2_key
                        if bin_compare not in sparse_matrix:
                            sparse_matrix[bin_compare] = [0, 0, 0, True]
                        elif sparse_matrix[bin_compare][3] == False:
                            continue
                        sparse_matrix[bin_compare][1] += 1
                    else:
                        break

    for bin_compare in sparse_matrix:
        chrom1, coord1, chrom2, coord2 = bin_compare.split(":")
        a, b, c, whitelist = sparse_matrix[bin_compare]
        if whitelist:
            denom = (b + c) / 2
            if denom > 0 and a > 0:
                statistic = a**2 / denom**2
                n = a + b + c
                lineHolder = (
                    "\t".join(
                        (
                            "0",
                            chrom1,
                            coord1,
                            "0",
                            "0",
                            chrom2,
                            coord2,
                            "1",
                            str(statistic),
                            str(a),
                            str(n),
                        )
                    )
                    + "\n"
                )
                out_lines.append(lineHolder)
            if len(out_lines) == 10000:
                with file_lock:
                    with open(output, "a") as out_file:
                        out_file.writelines(out_lines)
                        out_lines = []
    with file_lock:
        with open(output, "a") as out_file:
            out_file.writelines(out_lines)


def filter_short():
    short_file = "/home/bmt26/garg/fiberSeq/topology/bin_fire_mm10.short"
    output = "/home/bmt26/garg/fiberSeq/topology/bin_fire_mm10.>10.short"
    out_lines = []
    with open(short_file) as file:
        for line in file:
            cell = line.split("\t")
            if int(cell[-1]) > 10:
                out_lines.append(line)
            if len(out_lines) == 1000:
                with open(output, "a") as out_file:
                    out_file.writelines(out_lines)
                    out_lines = []
    with open(output, "a") as out_file:
        out_file.writelines(out_lines)


def powlaw(x, a, b, c):
    return a / numpy.power(x, b) + c


def max_fire_score_bedgraph():
    short_file = "/home/bmt26/garg/fiberSeq/topology/bin128i_fire_mm10.short"
    max_matrix = {}
    with open(short_file) as file:
        for line in file:
            cell = line.split("\t")
            if int(cell[10]) > 5:
                if abs(int(cell[2]) - int(cell[6])) != 128:
                    locus1 = cell[1] + "\t" + cell[2]
                    locus2 = cell[1] + "\t" + cell[6]
                    fire_score = float(cell[8])
                    if locus1 not in max_matrix:
                        max_matrix[locus1] = 0
                    if locus2 not in max_matrix:
                        max_matrix[locus2] = 0
                    if max_matrix[locus1] < fire_score:
                        max_matrix[locus1] = fire_score
                    if max_matrix[locus2] < fire_score:
                        max_matrix[locus2] = fire_score
    outlines = []
    for locus in max_matrix:
        fire_score = max_matrix[locus]
        chrom, center = locus.split("\t")
        i_center = int(center)
        start = str(i_center - 64)
        end = str(i_center + 64)
        out_line = chrom + "\t" + start + "\t" + end + "\t" + str(fire_score)
        outlines.append(out_line + "\n")
    with open("/home/bmt26/garg/fiberSeq/topology/bin128_fire_intergenic_max.bedgraph", "w") as file:
        file.writelines(outlines)


def mcfs_dist(dist_tupple):
    median_arr = []
    dist_arr = []
    x = []
    y = []
    for pair in dist_tupple:
        dist_arr.append(pair[0])
        median_arr.append(pair[1])
        if len(median_arr) == 2000:
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
    plt.savefig("/home/bmt26/garg/fiberSeq/topology/bin128d_fire_mm10.png")
    return (a, b, c)


def distance_correction():
    short_file = "/home/bmt26/garg/fiberSeq/topology/bin128f_fire_mm10.short"
    dist_tupple = []
    with open(short_file) as file:
        for line in file:
            cell = line.split("\t")
            # minimum n events
            if int(cell[9]) > 5:
                fire_dist = abs(int(cell[2]) - int(cell[6]))
                fire_score = float(cell[8])
                dist_tupple.append((fire_dist, fire_score))

    dist_tupple.sort(key=lambda tup: tup[0])
    mcfs_dist(dist_tupple)


def abcd_debugging():
    short_file = "/home/bmt26/garg/fiberSeq/topology/bin128f_fire_mm10.short"
    with open(short_file) as file:
        for line in file:
            cell = line.split("\t")
            locus1 = int(cell[2])
            locus2 = int(cell[6])
            if cell[1] == "chr3":
                if locus1 > 34759479 and locus1 < 34761938:
                    if locus2 > 34759479 and locus2 < 34761938:
                        print(locus1, locus2, cell[8])
    """
    json_dict = load_object(
        "/home/bmt26/garg/fiberSeq/topology/bin128d_fire_mm10_chr3.json"
    )
    for fiber in json_dict:
        elements_count = 0
        for locus in json_dict[fiber]:
            elements = json_dict[fiber][locus]
            pos = int(locus.split(":")[1])
            if pos > 34760450 and pos < 34760569:
                if int(elements) == 1:
                    elements_count += 1
            if pos > 34761480 and pos < 34761600:
                if int(elements) == 1:
                    elements_count += 1
        if elements_count > 0:
            print(elements_count, fiber)
    """


# abcd_debugging()

# distance_correction()

# Parts of program
binsToBed()
#bedtools intersect -sorted -a ${fire_peaks_only}.ep.bed -b ${fire_calls} -wo > ${intersect_output}
FIRE_intersect_to_matrix(input_bed=${intersect_output})
max_fire_score_bedgraph()
filter_short()
