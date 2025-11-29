bed_lines = []
gen_add = 0
survivors = 4
line_seq_holder = {}
import math
import matplotlib.pyplot as plt
import seaborn as sns
import os
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.manifold import TSNE


allowed_motifs = set(
    [
        "BACH2",
        "NFE2",
        "SOX2",
        "FOXC1",
        "RFX6",
        "FOXD3",
        "FOXA2",
        "SOX10",
        "MAFB",
        "RFX3",
        "ZIC3",
        "PO2F2",
        "MAFK",
        "ZBT17",
        "PO5F1",
        "NANOG",
        "LEF1",
        "PO3F1",
        "SOX3",
        "MAFG",
        "NF2L2",
        "SOX9",
        "BACH1",
        "PO2F1",
        "RFX2",
        "ZIC2",
        "RFX1",
        "REST",
        "NR5A2",
        "SOX4",
        "MAFF",
        "CTCF",
        "CTCFL",
        "SP4",
    ]
)

tsv_lines = []
# INSIDE_256_g117.out
with open("INSIDE_jun13.txt") as file:
    current_stat = 0
    gen_tracker = {}
    replicates_mean = {}
    itr = -1
    first_lines = True
    for line in file:
        line = line.strip()
        if "NEW_FILE" in line:
            gen_add += 1

        if "" not in line:
            continue

        if "Gen:" not in line:
            continue
        itr += 1
        surv_num = (itr % survivors) + 1
        cell = line.split("\t")
        cell[0] = "G" + str(int(cell[0].split("Gen:")[-1]) + gen_add)
        current_stat += float(cell[3])
        cell[1] += "_" + str(surv_num)
        if cell[1] not in line_seq_holder:
            line_seq_holder[cell[1]] = {}
        if cell[0] not in line_seq_holder[cell[1]]:
            line_seq_holder[cell[1]][cell[0]] = []
        line_seq_holder[cell[1]][cell[0]].append((cell[2], cell[3]))
# print(line_seq_holder)
with open("INSIDE_40_definitive/fimo.tsv") as file:
    for line in file:
        cell = line.split("\t")
        if cell[0] != "motif_id" and len(cell) > 5:
            i_line, gen, surv_rep, score = cell[2].split("_")
            # chr_n = ord(i_line[-1]) - ord("a") + 1
            chr_n = i_line.split("rep")[-1]
            chrom = f"chr{chr_n}"
            # print(gen)
            # print(cell)
            if float(cell[-3]) < 10**-6:
                # print(cell[0].split("_")[0])
                cell[0] = cell[0].split("_")[0] + "_" + score
                motif, score = cell[0].split("_")
                # motif = f"{motif}_{cell[5]}"

                b_line = "\t".join([chrom, cell[3], cell[4], cell[0], "1", cell[5]]) + "\n"
                line_name = i_line + "_" + surv_rep
                m_pos = int((int(cell[3]) + int(cell[4])) / 2)
                line_seq_holder[line_name][gen].append((motif, m_pos))
                # print(line_seq_holder[line_name][gen])
                # bed_lines.append(b_line)

ecdf_holder = {}
gen_skip = 1
samples_tfs = []
for line in line_seq_holder:
    for gen in line_seq_holder[line]:
        gen_mod = int(gen.replace("G", "")) % gen_skip
        # print(gen)
        if gen != "G30dhjfhhbf0":
            score = line_seq_holder[line][gen][0][1]
            seq = line_seq_holder[line][gen][0][0]
            motifs = line_seq_holder[line][gen][1:]
            motifs.sort(key=lambda x: x[1])
            if "_1" in line:
                if gen == "G94":
                    sample_tf = []
                    for m in motifs:
                        if m[0] in allowed_motifs:
                            sample_tf.append(m[0])
                    if len(sample_tf) != 0:
                        samples_tfs.append(sample_tf)

# red DO 1, 30, 60, 90, 120, 150
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import pairwise_distances


# Create a list of all unique TFs across all samples
all_tfs = sorted(set(tf for sample in samples_tfs for tf in sample))

# Create a binary matrix (samples x TFs), where 1 indicates presence and 0 indicates absence
sample_tf_matrix = np.array([[1 if tf in sample else 0 for tf in all_tfs] for sample in samples_tfs])

# Create a DataFrame for better handling with seaborn
df = pd.DataFrame(sample_tf_matrix, columns=all_tfs)

# Convert the DataFrame to a numpy array before calculating pairwise distances
df_array = df.values

# Calculate the Cosine distance matrix (1 - Cosine similarity)
cosine_dist_matrix = pairwise_distances(df_array.T, metric="jaccard")  # Calculate distance between TFs

# Plot the clustermap using the Cosine distance matrix with 'Reds' color palette
sns.clustermap(
    cosine_dist_matrix,
    figsize=(10, 8),
    cmap="Reds_r",
    method="ward",
    annot=False,
    row_cluster=True,
    col_cluster=True,
    xticklabels=all_tfs,
    yticklabels=all_tfs,
)

# Add title
# plt.title("Cosine Distance Clustermap of Transcription Factors")

# Show the plot
plt.show()
