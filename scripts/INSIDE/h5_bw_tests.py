import h5py
import os
import numpy
import math
import pybedtools

gen_obj = [['1','2','3','4'],['a','b','c','d'],['w','x','y','z']]


def score_predictions(gen_obj,score_regions,scored_parents=None):
    track = 1

    with h5py.File('predict.h5','r') as hf:
        syn_dnas = numpy.array(hf["chrom"][:])
        stat = numpy.ones(len(syn_dnas))
        preds = numpy.array(hf["preds"][:])
        for region in score_regions:
            bin_start = math.floor(region[0]/128)
            bin_end =  math.ceil(region[1]/128)
            stat *= preds[:,bin_start:bin_end,track].sum(axis=1)
    stat_itr = 0
    for line_i in range(len(gen_obj)):
        for seq_i in range(len(gen_obj[line_i])):
            seq = gen_obj[line_i][seq_i]
            gen_obj[line_i][seq_i] = (seq,stat[stat_itr])
            stat_itr += 1 
    if scored_parents == None:
        return gen_obj
    else:
        for line_i in range(len(scored_parents)):
            survivors_allowed = len(scored_parents[line_i]) 
            gen_obj[line_i] += scored_parents[line_i]
            gen_obj[line_i].sort(key=lambda x: x[1], reverse=True)
            print(gen_obj[line_i])
            gen_obj[line_i] = gen_obj[line_i][:survivors_allowed]
        return gen_obj

    

enh_len = 250
margin = 50
center_len = 900
seq_len = enh_len * 2 + margin * 4
print(seq_len)
center = "T"*center_len
half_loci_size = int( (131072 - (seq_len+center_len) )/ 2)
print(half_loci_size)
left = "T"* half_loci_size
right = "T"* half_loci_size

# tupples contain bed coords for regions to score
l_start = half_loci_size + margin
l_end = l_start + enh_len
r_start = l_end + margin + center_len + margin

score_regions = [(l_start,l_end),(r_start,r_start+enh_len)]
scored_parents = [[('1', .06251776218414307), ('2', 0.0304105281829834, b'l_1_s_2'), ('3', 0.5015015602111816, b'l_1_s_3'), ('4', 0.3099707365036011, b'l_1_s_4')], [('a', 0.412683898210525513, b'l_2_s_1'), ('b', 0.3120267391204834, b'l_2_s_2'), ('c', 0.9478468894958496, b'l_2_s_3'), ('d', 0.001802436076104641, b'l_2_s_4')], [('w', 0), ('x', 0.2977166175842285, b'l_3_s_2'), ('y', 0.317757415771484375, b'l_3_s_3'), ('z', 0.0016173720359802246, b'l_3_s_4')]]
print(score_predictions(gen_obj,score_regions,scored_parents))
