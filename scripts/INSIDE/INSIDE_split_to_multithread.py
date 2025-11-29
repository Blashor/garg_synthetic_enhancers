import json
import os

checkpoint = "latest-2.checkpoint"
j_path = "cluster140"
split = 4
# split has to be a clean multiple of the total


def read_checkpoint():
    with open(checkpoint) as json_file:
        checkpoint_json = json.load(json_file)
        scored_parents = checkpoint_json["scored_parents"]
        generation_start = checkpoint_json["current_gen"]
        fasta_info = checkpoint_json["fasta_info"]
        # print(fasta_info)
    return checkpoint_json


def read_setup_json(j_path):
    with open("setup_" + j_path + ".json") as json_file:
        setup_obj = json.load(json_file)
    return setup_obj


check = read_checkpoint()
setup = read_setup_json(j_path)
split_group_sizes = len(setup) / split

group_setup = {}
fasta_setup = {}
group_id = -1

for i, sample in enumerate(setup, start=0):
    if i % split_group_sizes == 0:
        group_id = i
        print(i, sample)
    if group_id not in group_setup:
        group_setup[group_id] = {}
    if group_id not in fasta_setup:
        fasta_setup[group_id] = {
            "current_gen": check["current_gen"],
            "scored_parents": [],
            "fasta_info": [],
        }

    group_setup[group_id][sample] = setup[sample]
    fasta_setup[group_id]["scored_parents"].append(check["scored_parents"][i])
    fasta_setup[group_id]["fasta_info"].append(check["fasta_info"][i])

# copy paste directories into INSIDE folder
# copy paste setups into json fold


os.system(f"mkdir -p {j_path}_split")
for group_id in group_setup:
    print(len(group_setup[group_id]))
    print(fasta_setup[group_id])
    with open(f"{j_path}_split/setup_{j_path}_{group_id}.json", "w") as file:
        file.write(json.dumps(group_setup[group_id]))
    os.system(f"mkdir -p {j_path}_split/{j_path}_{group_id}")
    with open(f"{j_path}_split/{j_path}_{group_id}/latest.checkpoint", "w") as file:
        file.write(json.dumps(fasta_setup[group_id]))
