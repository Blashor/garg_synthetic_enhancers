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


def output_scores(scored_parents, setup_obj, gen):
    line_itr = 0
    gen = str(gen + 1)
    for line_name in setup_obj:
        line_parents = scored_parents[line_itr]
        for parent in line_parents:
            print("Gen: " + gen, line_name, parent[0], str(parent[1]), sep="\t")
        line_itr += 1


def rollback(gen=150, checkpoint="latest.checkpoint"):
    # Restart from checkpoint
    gen_add = 0
    line_seq_holder = {}
    # INSIDE_256_g117.out
    with open("INSIDE_gen80_rep40.out") as file:
        current_stat = 0
        gen_tracker = {}
        replicates_mean = {}
        itr = 0
        first_lines = True
        for line in file:
            line = line.strip()
            if "NEW_FILE" in line:
                gen_add += 1

            if "" not in line:
                continue

            if "Gen:" not in line:
                continue

            cell = line.split("\t")
            cell[0] = "Gen: " + str(int(cell[0].split("Gen:")[-1]) + gen_add)
            current_stat += float(cell[3])
            if cell[1] not in line_seq_holder:
                line_seq_holder[cell[1]] = {}
            if cell[0] not in line_seq_holder[cell[1]]:
                line_seq_holder[cell[1]][cell[0]] = []
            line_seq_holder[cell[1]][cell[0]].append((cell[2], float(cell[3])))
    rollback_parents = []
    for sample in line_seq_holder:
        print(line_seq_holder[sample][f"Gen: {gen}"][1])
        rollback_parents.append(line_seq_holder[sample][f"Gen: {gen}"])
    # print(numpy.array(rollback_parents).shape)
    with open(checkpoint) as json_file:
        checkpoint_json = json.load(json_file)
    scored_parents = checkpoint_json["scored_parents"]
    generation_start = checkpoint_json["current_gen"]
    fasta_info = checkpoint_json["fasta_info"]
    print(generation_start)
    # print(numpy.array(scored_parents).shape)

    with open("rollback.checkpoint", "w") as file:
        file.write(
            json.dumps(
                {
                    "scored_parents": rollback_parents,
                    "current_gen": gen,
                    "fasta_info": fasta_info,
                }
            )
        )

    sys.stdout.flush()


def read_setup_json(j_path):
    with open("setup_" + j_path + ".json") as json_file:
        setup_obj = json.load(json_file)
    return setup_obj


setup_obj = read_setup_json("cluster140")

rollback()

# setup_obj = setup_obj_helper()
# main(setup_obj, generations=100)
