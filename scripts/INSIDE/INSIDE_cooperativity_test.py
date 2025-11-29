import random
import os
import sys
import numpy
import h5py
import math
import pybedtools
import json
import copy


def snp_mutator(seq, num_mutations, child_sequences):
    # print("Seq Len:", len(seq))
    child_seqs = []
    seq_len = len(seq)
    for child_i in range(child_sequences):
        child_seq = list(seq)
        positions = random.sample(range(seq_len), num_mutations)
        for position in positions:
            bases = ["A", "C", "G", "T"]
            bases.remove(seq[position])
            new_base = random.choice(bases)
            child_seq[position] = new_base
        child_seqs.append("".join(child_seq))
    return child_seqs


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
                if seq_itr > 1:
                    seq_left = snp_mutator(seq_left, int(0.75 * len(seq_left)), 1)[0]
                seq = seq_left + center + seq_right
            fasta_lines.append(left + seq + right + "\n")
            bed_lines.append(seq_id + "\t0\t131072\t\n")
    with open(working_dir + "INSIDE_seqs_test.fa", "w") as file:
        file.writelines(fasta_lines)
    os.system("samtools faidx " + working_dir + "INSIDE_seqs_test.fa")
    with open(working_dir + "INSIDE_seqs_test.bed", "w") as file:
        file.writelines(bed_lines)


def coop_test(checkpoint):
    with open(checkpoint) as json_file:
        checkpoint_json = json.load(json_file)
    scored_parents = checkpoint_json["scored_parents"]
    generation_start = checkpoint_json["current_gen"]
    fasta_info = checkpoint_json["fasta_info"]
    gen_obj = []
    for line in scored_parents:
        new_line = [line[0][0], line[1][0], line[0][0], line[1][0]]
        gen_obj.append(new_line)

    # Get scores for current gen children
    seqs_to_bedANDfasta(gen_obj, fasta_info)
    run_prediction()


def run_prediction():
    fasta = working_dir + "INSIDE_seqs_test.fa"

    wig_table = "../data/" + model_name + "/targets.txt"

    trained_model = "../models/" + model_name
    model_params = trained_model + "/params.json"

    test_output = working_dir + model_name + "_test"

    bed = working_dir + "INSIDE_seqs_test.bed"

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


def get_predicts(checkpoint):
    with open(checkpoint) as json_file:
        checkpoint_json = json.load(json_file)
    scored_parents = checkpoint_json["scored_parents"]
    generation_start = checkpoint_json["current_gen"]
    fasta_info = checkpoint_json["fasta_info"]
    track = 1
    stat_obj = {"left_delta": [], "right_delta": []}
    with h5py.File(working_dir + model_name + "_test/" + "predict.h5", "r") as hf:
        syn_dnas = numpy.array(hf["chrom"][:])
        stat = numpy.ones(len(syn_dnas))
        preds = numpy.array(hf["preds"][:])
    print(preds.shape)
    fa_itr = 0
    for line_itr in range(len(scored_parents)):
        score_regions = fasta_info[line_itr][-1]
        sur_stat_obj = {}
        for survivor in scored_parents[line_itr]:
            stat_arr = []
            for region in score_regions:
                bin_start = math.floor(region[0] / 128)
                bin_end = math.ceil(region[1] / 128)
                stat_arr.append(numpy.mean(preds[fa_itr, bin_start:bin_end, track]))
            print(syn_dnas[fa_itr], stat_arr, fa_itr % 4)
            if fa_itr % 4 == 0:
                sur_stat_obj["left1"] = stat_arr[0]
                sur_stat_obj["right1"] = stat_arr[1]
            if fa_itr % 4 == 1:
                sur_stat_obj["left2"] = stat_arr[0]
                sur_stat_obj["right2"] = stat_arr[1]
            if fa_itr % 4 == 2:
                if sur_stat_obj["left1"] > 5 and sur_stat_obj["right1"] > 5:
                    sur_stat_obj["left1"] -= stat_arr[0]
                    sur_stat_obj["right1"] -= stat_arr[1]
                else:
                    sur_stat_obj["left1"] = 0
                    sur_stat_obj["right1"] = 0
            if fa_itr % 4 == 3:
                if sur_stat_obj["left2"] > 5 and sur_stat_obj["right2"] > 5:
                    sur_stat_obj["left2"] -= stat_arr[0]
                    sur_stat_obj["right2"] -= stat_arr[1]
                else:
                    sur_stat_obj["left2"] = 0
                    sur_stat_obj["right2"] = 0
                stat_obj["left_delta"].append(sur_stat_obj["left1"])
                stat_obj["left_delta"].append(sur_stat_obj["left2"])
                stat_obj["right_delta"].append(sur_stat_obj["right1"])
                stat_obj["right_delta"].append(sur_stat_obj["right2"])
            fa_itr += 1
    stat_obj["left_delta"] = numpy.array(stat_obj["left_delta"])
    stat_obj["right_delta"] = numpy.array(stat_obj["right_delta"])

    stat_obj["left_delta"] = stat_obj["left_delta"][stat_obj["left_delta"].nonzero()]
    stat_obj["right_delta"] = stat_obj["right_delta"][stat_obj["right_delta"].nonzero()]
    print(stat_obj["left_delta"])
    print(stat_obj["right_delta"])
    print(numpy.median(stat_obj["left_delta"]), numpy.mean(stat_obj["left_delta"]))
    print(numpy.median(stat_obj["right_delta"]), numpy.mean(stat_obj["right_delta"]))


working_dir = "../INSIDE/INSIDE_mut/"
model_name = "ce_concord_v4_3"

coop_test(working_dir + "latest.checkpoint")
get_predicts(
    working_dir + "latest.checkpoint",
)

"""
get_predicts(
    working_dir + "latest.checkpoint",
)
"""
