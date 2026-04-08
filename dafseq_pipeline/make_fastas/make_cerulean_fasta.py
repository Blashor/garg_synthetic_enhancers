# this will read dna files
# create fastas for each corresponding demux from dna file
# align

from snapgene_reader import snapgene_file_to_dict, snapgene_file_to_seqrecord
import os

("chr1", 2448, 2964), ("chr1", 2964, 3482)

dna_dir = "/Users/blake/Documents/gargLab/INSIDE_atac/demux_files/DNA_files"
snaps = os.listdir(dna_dir)

f1 = "G1."
f2 = "G60"
lines = []
for i, snap in enumerate(snaps):
    fa_name = snap.replace(".dna", "")
    seq = snapgene_file_to_dict(f"{dna_dir}/{snap}")["seq"]
    lines.append(f">{fa_name}\n{seq}\n")
with open("/Users/blake/Documents/gargLab/INSIDE_atac/demux_files/r6_insert.fa", "w") as file:
    file.writelines(lines)
