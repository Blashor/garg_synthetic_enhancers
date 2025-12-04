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
def chrome_parser(chrom_str):
    chrom, coords = chrom_str.split(":")
    start, end = map(int, coords.split("-"))
    return chrom, start, end


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
                seq_left, seq_right = (seq[0], seq[1])
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
        score_regions = fasta_info[line_i][-1]
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


def region_setup(l, r):
    chrom, start_l, end_l = chrome_parser(l)
    chrom, start_r, end_r = chrome_parser(r)
    center_len = start_r - end_l
    enh_l = end_l - start_l
    enh_r = end_r - start_r

    seq_len = enh_l + enh_r  # __left__mEm_center_mEm__right__
    half_loci_start = int((131072 - (seq_len + center_len)) / 2)

    l_start = half_loci_start
    l_end = l_start + enh_l
    r_start = l_end + center_len
    score_regions = [(l_start, l_end), (r_start, r_start + enh_r)]

    left_side = (chrom, start_l - score_regions[0][0], start_l)  # left region size
    left_enh = (chrom, start_l, end_l)
    center_region = (chrom, end_l, start_r)
    right_enh = (chrom, start_r, end_r)
    right_side = (chrom, end_r, end_r + (131072 - score_regions[1][1]))
    regions = (left_side, right_side, center_region, left_enh, right_enh)
    seq_obj = []
    for region in regions:
        seq_obj.append(str(pybedtools.BedTool.seq(region, genome)))

    seq_obj.append((enh_l, enh_r))  # sizes
    seq_obj.append(score_regions)
    return seq_obj


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
            fasta_info.append(region_setup(line_info["enh_l"], line_info["enh_r"]))

            # Prediction Calculations
            total_predictions += line_info["survivors"]
            total_predictions += (
                generations * line_info["children_per"] * line_info["survivors"] * len(line_info["mut_rates"])
            )

            # Create Random Sequences
            enh_l = fasta_info[-1][3]
            enh_r = fasta_info[-1][4]
            survs_arr = []
            for s in range(line_info["survivors"]):
                survs_arr.append([snp_mutator(enh_l, 0, 1)[0], snp_mutator(enh_r, 0, 1)[0]])

            gen_zero_obj.append(survs_arr)

        print("Total number of Predictions:", total_predictions)
        print(setup_obj)
        print(gen_zero_obj)
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
            # seq_len = fasta_info[line_itr][3]

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
                            child_seq_l = snp_mutator(parent[0], mut_obj["SNPs"], 1)[0]
                            child_seq_r = snp_mutator(parent[1], mut_obj["SNPs"], 1)[0]
                            child_seq = (child_seq_l, child_seq_r)
                            # if mut_obj["dups"] > 0:
                            #    child_seq = segment_copier(child_seq, 1, mut_obj["dups"])[0]
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


# sys.argv = ["python", "supersInGenome"]
setup_obj = read_setup_json(sys.argv[1])

working_dir = "../INSIDE/" + sys.argv[1] + "/"
model_name = "ce_concord_v4_3"
model_name = "cluster_mar_28"
os.system("mkdir -p " + working_dir)

genome = pybedtools.BedTool("/gpfs/gibbs/project/garg/shared/catalogs/gencode_fastas/GRCm38.primary_assembly.genome.fa")
# setup_obj = setup_obj_helper()
# main(setup_obj, generations=100)
try:
    main(setup_obj, checkpoint=working_dir + "latest.checkpoint", generations=5000)
except FileNotFoundError:
    main(setup_obj, generations=5000)
