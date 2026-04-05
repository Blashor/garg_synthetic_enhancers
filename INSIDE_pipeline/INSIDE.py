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


def mutation_bag(chance, num_times):
    mut_obj = {"SNPs": 0, "dups": 0}
    for i in range(num_times):
        res = random.random() < chance
        if res == True:
            mut_obj["dups"] += 1
        else:
            mut_obj["SNPs"] += 1
    return mut_obj


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


def segment_copier(seq, child_sequences, duplication_mutations):
    seq_size = len(seq)
    child_seqs = []
    for child_i in range(child_sequences):
        # selects a segment to copy
        segment_size = random.randint(4, 16)
        upper_bound = seq_size - (segment_size + 1)
        copy_site = random.randint(0, upper_bound)
        copy_seq = seq[copy_site : copy_site + segment_size]

        # num of times segment is duplicated
        copy_num = duplication_mutations
        child_seq = seq
        for copy_i in range(copy_num):
            insert_site = random.randint(0, upper_bound)
            child_seq = child_seq[0:insert_site] + copy_seq + child_seq[insert_site + segment_size :]
        child_seqs.append(child_seq)
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
            gen_obj[line_i] += scored_parents[line_i]
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
    if checkpoint == None:
        total_predictions = 0
        for line_name in setup_obj:
            # Setup
            line_info = setup_obj[line_name]
            fasta_info.append(region_setup(line_info["enh_len"], line_info["margin"], line_info["center_len"]))

            # Prediction Calculations
            total_predictions += line_info["survivors"]
            total_predictions += (
                generations * line_info["children_per"] * line_info["survivors"] * len(line_info["mut_rates"])
            )

            # Create Random Sequences
            seq_len = fasta_info[-1][3]
            gen0_parents = snp_mutator("T" * seq_len, int(seq_len * 0.75), line_info["survivors"])
            gen_zero_obj.append(gen0_parents)

        print("Total number of Predictions:", total_predictions)
        print(setup_obj)

        # Get inital scores for generation 0
        seqs_to_bedANDfasta(gen_zero_obj, fasta_info)
        run_prediction()
        scored_parents = score_predictions(gen_zero_obj, fasta_info)
        generation_start = 0
    else:
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
                for mut in line_info["mut_rates"]:
                    child_seqs = []
                    if mut >= 1:
                        num_muts = int(math.floor(mut))
                        dup_chance = mut - math.floor(mut)
                        mut_obj = mutation_bag(dup_chance, num_muts)
                        for child_i in range(line_info["children_per"]):
                            child_seq = snp_mutator(parent, mut_obj["SNPs"], 1)[0]
                            if mut_obj["dups"] > 0:
                                child_seq = segment_copier(child_seq, 1, mut_obj["dups"])[0]
                            child_seqs.append(child_seq)
                    line_pool += child_seqs
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
    with open("../INSIDE/json/setup_" + j_path + ".json") as json_file:
        setup_obj = json.load(json_file)
    return setup_obj


setup_obj = read_setup_json(sys.argv[1])

working_dir = "../INSIDE/" + sys.argv[1] + "/"
model_name = "ce_concord_v4_3"
model_name = "cluster_mar_28"
os.system("mkdir -p " + working_dir)

# setup_obj = setup_obj_helper()
# main(setup_obj, generations=100)
try:
    main(setup_obj, checkpoint=working_dir + "latest.checkpoint", generations=150)
except FileNotFoundError:
    main(setup_obj, generations=5000)
