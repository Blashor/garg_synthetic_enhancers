import json
import sys


def prep_object(prefix):
    lines = []
    # opens the bed file where each line has: SE, CE, Fiber ID, FDR
    # print(prefix)
    with open(prefix + "/Cov.bed", "r") as file:
        lines = file.read().split("\n")
    for line in lines:
        fiber_id = line[: line.rfind("\t")]
        # checks FDR of each line, returns a dict in format of se_fire_score[super_id][enh_id][fiber_id] = {"best_score":,"length":}
        fiber_to_se(line)

    # This moves it into a format that gets total local fiber count so a data matrix can be made
    matrix = []
    for se_id, se in se_fire_score.items():
        if se_id != "":
            se_mat = {"seId": se_id, "enhs": []}
            for enh_id, enh in se.items():
                if enh_id != "fiberList":
                    enh_range = enh_id.split(":")[1].split("-")
                    enh_size = int(enh_range[1]) - int(enh_range[0])
                    enh_column = {"enhId": enh_id, "fibers": []}
                    for fiber in se["fiberList"]:
                        # goes through each fiber in se_fire_score, using "fiberList" ensures that index is shared by the same fiber across enhancers
                        try:
                            fire_score = float(enh[fiber]["best_score"])
                            if fire_score < 0.1:
                                enh_column["fibers"].append(fire_score)
                            elif enh[fiber]["length"] == enh_size:
                                enh_column["fibers"].append(fire_score)
                            else:
                                enh_column["fibers"].append(None)
                        except KeyError:
                            enh_column["fibers"].append(None)
                    se_mat["enhs"].append(enh_column)
            matrix.append(se_mat)
    with open(prefix + "_obj.json", "w") as file:
        file.write(json.dumps(matrix))
    return matrix


def fiber_to_se(fiber):
    cells = fiber.split("\t")
    if len(cells) == 9:
        enh_id = cells[0] + ":" + cells[1] + "-" + cells[2]
        super_id = cells[3] + ":" + cells[4] + "-" + cells[5]
        fiber_id = cells[6]
        fire_score = cells[7]
        feature_length = int(cells[8])
        # print(cells)
        if super_id not in se_fire_score:
            se_fire_score[super_id] = {"fiberList": []}
        if enh_id not in se_fire_score[super_id]:
            se_fire_score[super_id][enh_id] = {}
        if fiber_id not in se_fire_score[super_id]["fiberList"]:
            se_fire_score[super_id]["fiberList"].append(fiber_id)
        # this is to deal with nucleosomes + FIRE elements both overlapping underneath a peak
        if fiber_id in se_fire_score[super_id][enh_id]:
            if fire_score < se_fire_score[super_id][enh_id][fiber_id]["best_score"]:
                se_fire_score[super_id][enh_id][fiber_id]["best_score"] = fire_score
            se_fire_score[super_id][enh_id][fiber_id]["length"] += feature_length
        else:
            se_fire_score[super_id][enh_id][fiber_id] = {
                "best_score": fire_score,
                "length": feature_length,
            }


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


def main(prefix):
    print(prefix)
    global lengthCutOff
    global se_fire_score
    se_fire_score = {}
    lengthCutOff = 0
    covArrs = []
    mat = prep_object(prefix)
    # mat = load_object()
    # print(len(mat))
    # covArrs.append(get_correlation(mat))
    # print(scipy.stats.ks_2samp(covArrs[0], covArrs[1]))


# main()


lengths = [
    0,
    10,
    100,
    500,
    1000,
    5000,
    10000,
]
