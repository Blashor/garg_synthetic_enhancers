import random
import os
import sys
import numpy
import h5py
import math
import pybedtools
import json
import copy


def setup_obj_helper(file_name, new_enh_lens=[], new_center_lens=[], mut_rates=[], margins=[]):
    setup_obj = {}
    example_line = {
        "survivors": 4,
        "children_per": 500,
        "mut_rates": [2.0],  # new setup of number of mutation, percent that are dups
        "enh_len": 256,
        "margin": 0,
        "center_len": 1024,
    }
    replicates = list(map(lambda x: "rep" + str(x + 1), range(3)))
    #
    #
    # For mut rate adjustment based on enh size

    for rep in replicates:
        for e_len in new_enh_lens:
            max_e_len = max(new_enh_lens)
            max_children_per = example_line["children_per"]
            line_name = "line" + str(e_len) + rep
            new_line = copy.deepcopy(example_line)
            new_line["enh_len"] = e_len
            # new_line["children_per"] = round(max_children_per * (e_len / max_e_len))
            setup_obj[line_name] = new_line
        for c_len in new_center_lens:
            line_name = "line_center_" + str(c_len) + rep
            new_line = copy.deepcopy(example_line)
            if c_len == 0:
                new_line["center_len"] = 0
            new_line["center_len"] = c_len
            setup_obj[line_name] = new_line
        for m_rate in mut_rates:
            line_name = "line_mut_" + str(m_rate) + rep
            new_line = copy.deepcopy(example_line)
            new_line["mut_rates"] = [m_rate]
            setup_obj[line_name] = new_line
        for margin in margins:
            line_name = "line_margin_" + str(margin) + rep
            new_line = copy.deepcopy(example_line)
            new_line["margin"] = margin
            setup_obj[line_name] = new_line
    with open(file_name, "w") as file:
        file.write(json.dumps(setup_obj))
    for i in setup_obj:
        print(i, setup_obj[i])
    return setup_obj


# new_enh_lens = [2000, 1000, 500, 450, 400, 350, 300, 250, 200, 150, 100, 50, 25]
# new_center_lens = [0, 10, 100, 350, 1000, 2000, 3750, 7500, 12500, 15000]
# margins = [0, 25, 50, 75, 100, 125, 150, 175, 200, 225]


setup_obj_helper("json/setup_mutations.json", mut_rates=[1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0])
setup_obj_helper("json/setup_centers.json", new_center_lens=[128, 256, 512, 1024, 2048, 4096])
setup_obj_helper("json/setup_margins.json", margins=[0, 64, 128, 256, 512])
setup_obj_helper("json/setup_lens.json", new_enh_lens=[64, 128, 256, 512, 1024, 2048])
# one_ce_obj = setup_obj_helper("setup_one500Kids.json", new_enh_lens=new_enh_lens)
# e_len_obj = setup_obj_helper("setup_enh.json", new_enh_lens=new_enh_lens)
# margin_obj = setup_obj_helper("setup_margin.json", margins=margins)
# center_obj = setup_obj_helper("setup_center.json", new_center_lens=new_center_lens)
