# lets make something to read ce_ranks file and get images of top ones
import os
import sys
import json
import itertools


def run_sat(sat_name, sat_len):
    fasta = working_dir + "INSIDE_seqs.fa"
    #
    wig_table = "../data/" + model_name + "/targets.txt"
    trained_model = "../models/" + model_name
    model_params = trained_model + "/params.json"
    #
    sat_output = working_dir + "sat_output"

    sat_bed = working_dir + "sat_beds/" + sat_name + ".bed"
    sat_cmd = f"python basenji_sat_bed.py -f {fasta} -l {sat_len} -o {sat_output}/{sat_name} --rc -t {wig_table} {model_params} {trained_model}/model_best.h5 {sat_bed}"

    sat_plt_cmd = f"python basenji_sat_plot.py --png -l {sat_len} -o {sat_output}/{sat_name}/plots -t {wig_table} {sat_output}/{sat_name}/scores.h5"

    print(sat_cmd)
    os.system(sat_cmd)
    sys.stdout.flush()
    print(sat_plt_cmd)
    sys.stdout.flush()
    os.system(sat_plt_cmd)


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


def parse_parents_filter(p_file, seqs_to_test):
    line_ids = []
    parent_seqs = []
    with open(p_file) as file:
        for line in file:
            gen, line_name, seq, score = line.strip().split("\t")
            l_id = gen + "\t" + line_name
            if l_id in seqs_to_test:
                # tests first replicate only
                seqs_to_test.remove(l_id)
                l_id = l_id.replace("Gen: ", "g").replace("\t", "_")
                line_ids.append(l_id)
                parent_seqs.append(seq)
    return (line_ids, [parent_seqs])


def seqs_to_bedsANDfasta(next_generation_obj, line_ids, fasta_info):
    line_itr = 0
    fasta_lines = []

    for line in next_generation_obj:
        left, right, center = fasta_info[line_itr][:3]
        line_itr += 1
        seq_itr = 0
        for seq, seq_id in zip(line, line_ids):
            seq_itr += 1

            if center != None:
                seq_left, seq_right = seq[: len(seq) // 2], seq[len(seq) // 2 :]
                # seq = seq_left + center + seq_right
                skew_left = len(center) / 2 + len(seq_left) / 2
                skew_right = len(center) / 2 + len(seq_right) / 2
                #
                l_skew_left = "T" * int(len(left) + skew_left)
                l_skew_right = "T" * int(len(right) - skew_left)
                r_skew_left = "T" * int(len(left) - skew_right)
                r_skew_right = "T" * int(len(right) + skew_right)
                # print(len(l_skew_left), len(seq_left), len(center), len(seq_right), len(l_skew_right))
                # print(len(r_skew_left), len(seq_left), len(center), len(seq_right), len(r_skew_right))
                left_seq = l_skew_left + seq_left + center + seq_right + l_skew_right
                right_seq = r_skew_left + seq_left + center + seq_right + r_skew_right
            bed_lines = []
            fasta_lines.append(">" + seq_id + "_left\n")
            fasta_lines.append(left_seq + "\n")
            bed_lines.append(seq_id + "_left\t0\t131072\n")
            #
            fasta_lines.append(">" + seq_id + "_right\n")
            fasta_lines.append(right_seq + "\n")
            bed_lines.append(seq_id + "_right\t0\t131072\n")
            with open(working_dir + "sat_beds/" + seq_id + ".bed", "w") as file:
                file.writelines(bed_lines)
    with open(working_dir + "INSIDE_seqs.fa", "w") as file:
        file.writelines(fasta_lines)
    os.system("samtools faidx " + working_dir + "INSIDE_seqs.fa")


def main(setup_obj, parent_file, seqs_to_test):
    # Array of arrays, where each first order array is a line, and each seq in each array is a parent
    fasta_info = []
    # Initialize parents_obj
    for line_name in setup_obj:
        # Setup #doing this once works since I don't change anything between lines
        line_info = setup_obj[line_name]
        fasta_info.append(region_setup(line_info["enh_len"], line_info["margin"], line_info["center_len"]))
    # print(fasta_info)
    line_ids, parent_seqs = parse_parents_filter(parent_file, seqs_to_test)

    seqs_to_bedsANDfasta(parent_seqs, line_ids, fasta_info)
    sat_len = 512
    for sat_name in line_ids:
        run_sat(sat_name, sat_len)


def read_setup_json(j_path):
    with open("../INSIDE/json/setup_" + j_path + ".json") as json_file:
        setup_obj = json.load(json_file)
    return setup_obj


sys.argv.append("cluster140")
setup_obj = read_setup_json(sys.argv[1])

# model_name = "ce_concord_v4_3"
model_name = "cluster_mar_28"

working_dir = "../INSIDE/" + sys.argv[1] + "_sat_mutagenesis_" + model_name + "/"
parents_file = "/home/bmt26/scripts/INSIDE/INSIDE_jun13.txt"
# parents_file = "INSIDE_jun13.txt"
os.system("mkdir -p " + working_dir)
os.system("mkdir -p " + working_dir + "sat_beds")
os.system("mkdir -p " + working_dir + "sat_output")
samples = ["line256rep6", "line256rep28", "line256rep31"]
gens = list(map(lambda x: f"Gen: {x}", [1, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300]))
#
seqs_to_test = set(map(lambda x: "\t".join(x), list(itertools.product(gens, samples))))


#

#
#
main(setup_obj, parents_file, seqs_to_test)
