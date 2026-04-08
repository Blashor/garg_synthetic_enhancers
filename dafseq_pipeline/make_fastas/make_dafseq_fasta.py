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
for i, snap in enumerate(snaps):
    if f1 in snap:
        f1seq = snapgene_file_to_dict(f"{dna_dir}/{snap}")["seq"]
        print(i, snap)
    if f2 in snap:
        f2seq = snapgene_file_to_dict(f"{dna_dir}/{snap}")["seq"]
        print(i, snap)

f1seq = list(f1seq)
f1seq[2964:3482] = f2seq[2964:3482]
seq = "".join(f1seq)

f1 = f1.replace(".", "")
f2 = f2.replace(".", "")
with open(f"/Users/blake/Documents/gargLab/INSIDE_atac/demux_files/fake_fastas/{f1}_{f2}.fasta", "w") as file:
    file.write(f">chr1\n{seq}")
