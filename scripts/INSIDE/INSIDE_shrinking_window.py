import random
import os
import sys
import numpy
import h5py
import math
import pybedtools
import json
import copy

###
# INSIDE.py
#
# in silico directed evolution of sequence via basneji selection
#
# 1. Generates random sequences (including synthetic fastas + .bed to access)
# 2. Call basenji_predict_bed.py to score sequences
# 3. Select best children -> Repeat
###


def seqs_to_bedANDfasta(next_generation_obj, fasta_info):
    line_itr = 0
    fasta_lines = []
    bed_lines = []
    for line in next_generation_obj:
        left, right, center = fasta_info[line_itr][:3]
        line_itr += 1
        seq_itr = 0
        for seq in line:
            seq_itr += 1
            seq_id = "l_" + str(line_itr) + "_s_" + str(seq_itr)
            fasta_lines.append(">" + seq_id + "\n")
            if center != None:
                seq_left, seq_right = seq[: len(seq) // 2], seq[len(seq) // 2 :]
                seq = seq_left + center + seq_right
            fasta_lines.append(left + seq + right + "\n")
            bed_lines.append(seq_id + "\t0\t131072\t\n")
    with open(working_dir + "INSIDE_seqs.fa", "w") as file:
        file.writelines(fasta_lines)
    os.system("samtools faidx " + working_dir + "INSIDE_seqs.fa")
    with open(working_dir + "INSIDE_seqs.bed", "w") as file:
        file.writelines(bed_lines)


def run_prediction():
    fasta = working_dir + "INSIDE_seqs.fa"

    wig_table = "../data/" + model_name + "/targets.txt"

    trained_model = "../models/" + model_name
    model_params = trained_model + "/params.json"

    test_output = working_dir + model_name

    bed = working_dir + "INSIDE_seqs.bed"

    pred_cmd = (
        "python basenji_predict_bed.py -f "
        + fasta
        + " -o "
        + test_output
        + " --rc -t "
        + wig_table
        + " "
        + model_params
        + " "
        + trained_model
        + "/model_best.h5 "
        + bed
        + " > /dev/null"
    )

    # print(pred_cmd)
    os.system(pred_cmd)
    sys.stdout.flush()


def shrinking_window(seq):
    seq_left, seq_right = seq[: len(seq) // 2], seq[len(seq) // 2 :]
    new_seqs = []

    for small_seq, no_change, reverse, flip_seqs in [
        (seq_left, seq_right, False, False),
        (seq_left, seq_right, True, False),
        (seq_right, seq_left, False, True),
        (seq_right, seq_left, True, True),
    ]:
        shrink_seq = ""
        ran_into_T = True
        if reverse == True:
            small_seq = small_seq[::-1]
        for i, base in enumerate(small_seq):
            if ran_into_T == True:
                if base == "T":
                    # print(i, base)
                    shrink_seq += base
                else:
                    # print(i, "T")
                    shrink_seq += "T"
                    ran_into_T = False
            elif ran_into_T == False:
                # print(i, base)
                shrink_seq += base
        if reverse == True:
            shrink_seq = shrink_seq[::-1]
        if flip_seqs == False:
            # print(shrink_seq, no_change)
            new_seqs.append(shrink_seq + no_change)
        else:
            # print(no_change, shrink_seq)
            new_seqs.append(no_change + shrink_seq)
    # print(new_seqs)
    return new_seqs


def score_predictions(gen_obj, fasta_info, scored_parents=None):
    track = 1
    stat_itr = 0
    # Read basenji prediction file
    with h5py.File(working_dir + model_name + "/predict.h5", "r") as hf:
        syn_dnas = numpy.array(hf["chrom"][:])
        preds = numpy.array(hf["preds"][:])

    # Score regions
    for line_i in range(len(gen_obj)):
        score_regions = fasta_info[line_i][4]
        for seq_i in range(len(gen_obj[line_i])):
            seq = gen_obj[line_i][seq_i]
            stat = 1
            for region in score_regions:
                bin_start = math.floor(region[0] / 128)
                bin_end = math.ceil(region[1] / 128)
                stat *= numpy.mean(preds[stat_itr, bin_start:bin_end, track])
            gen_obj[line_i][seq_i] = (seq, stat)
            stat_itr += 1
    if scored_parents == None:
        return gen_obj
    else:
        for line_i in range(len(scored_parents)):
            survivors_allowed = len(scored_parents[line_i])
            # gen_obj[line_i] += scored_parents[line_i]
            gen_obj[line_i].sort(key=lambda x: x[1], reverse=True)
            gen_obj[line_i] = gen_obj[line_i][:survivors_allowed]
        return gen_obj


def output_scores(scored_parents, setup_obj, gen):
    line_itr = 0
    gen = str(gen + 1)
    for line_name in setup_obj:
        line_parents = scored_parents[line_itr]
        for parent in line_parents:
            print("Gen: " + gen, line_name, parent[0], str(parent[1]), sep="\t")
        line_itr += 1


def region_setup(enh_len, margin, center_len):
    seq_len = enh_len * 2 + margin * 4  # __left__mEm_center_mEm__right__
    half_loci_size = int((131072 - (seq_len + center_len)) / 2)
    left = "T" * half_loci_size
    right = "T" * half_loci_size
    center = "T" * center_len

    l_start = half_loci_size + margin
    l_end = l_start + enh_len
    r_start = l_end + margin + center_len + margin
    score_regions = [(l_start, l_end), (r_start, r_start + enh_len)]

    return (left, right, center, seq_len, score_regions)


def main(setup_obj, generations=50, checkpoint=None):
    # Array of arrays, where each first order array is a line, and each seq in each array is a parent
    gen_zero_obj = []
    fasta_info = []
    # Initialize parents_obj
    # Restart from checkpoint
    with open(checkpoint) as json_file:
        checkpoint_json = json.load(json_file)
    scored_parents = checkpoint_json["scored_parents"]
    generation_start = checkpoint_json["current_gen"]
    fasta_info = checkpoint_json["fasta_info"]

    print(scored_parents)
    sys.stdout.flush()
    # Run directed evolution
    for generation in range(generation_start, generations):
        next_generation_obj = []
        line_itr = 0
        # Create all line pools for a generation
        for line_name in setup_obj:
            # Setup
            line_info = setup_obj[line_name]
            parent_line = scored_parents[line_itr]
            seq_len = fasta_info[line_itr][3]

            # Survivors become parents => mutated children generated
            line_pool = []
            for parent, score in parent_line:
                child_seqs = []
                shrinking_children = shrinking_window(parent)
                line_pool += shrinking_children
            next_generation_obj.append(line_pool)

            line_itr += 1

        # Get scores for current gen children
        seqs_to_bedANDfasta(next_generation_obj, fasta_info)
        run_prediction()
        scored_parents = score_predictions(next_generation_obj, fasta_info, scored_parents)

        with open(working_dir + "latest.checkpoint", "w") as file:
            file.write(
                json.dumps(
                    {
                        "scored_parents": scored_parents,
                        "current_gen": generation,
                        "fasta_info": fasta_info,
                    }
                )
            )

        output_scores(scored_parents, setup_obj, generation)
        sys.stdout.flush()


def read_setup_json(j_path):
    with open("../INSIDE_shrinking_window_4_3/json/setup_" + j_path + ".json") as json_file:
        setup_obj = json.load(json_file)
    return setup_obj


setup_obj = read_setup_json(sys.argv[1])

working_dir = "../INSIDE_shrinking_window_4_3/" + sys.argv[1] + "/"
model_name = "ce_concord_v4_3"
# model_name = "cluster_mar_28"
os.system("mkdir -p " + working_dir)

# setup_obj = setup_obj_helper()
# main(setup_obj, generations=100)
try:
    main(setup_obj, checkpoint=working_dir + "latest.checkpoint", generations=5000)
except FileNotFoundError:
    main(setup_obj, generations=5000)
