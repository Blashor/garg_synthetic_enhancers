import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import os
import seaborn as sns
import math
import matplotlib.cm as cm
import random

#
# networkx style visual of motif occurences
#


def rearrange_adjacent(sorted_arr):
    n = len(sorted_arr)
    result = [None] * n
    mid = n // 2
    result[::2] = sorted_arr[:mid]
    result[1::2] = sorted_arr[mid:]
    return result


def mHstr(motif_holder):
    m_h = []
    for motif in motif_holder:
        motif_id = f"{motif}_{motif_holder[motif]}"
        m_h.append(motif_id)
    return m_h


# Create graph and add edges with weights

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
        "NR5A2",
        "SOX4",
        "MAFF",
        "CTCF",
        "CTCFL",
        "None",
    ]
)

allowed_motifs = set(["BACH2", "SOX2", "ZIC3", "PO2F2", "PO5F1", "NANOG", "LEF1", "PO3F1", "PO2F1", "CTCF", "None"])


def net_x(graphing_obj, node_scores, motif_sizes, generations, motif_appears_in_replicates):
    unique_motifs = set()
    G = nx.Graph()

    # Add edges based on the graphing object
    for new_motif in graphing_obj:
        for old_motif in graphing_obj[new_motif]:
            if node_scores[new_motif] < 1500 and node_scores[old_motif] < 1500:
                # graphing_obj[new_motif][old_motif] > 40 or
                # if graphing_obj[new_motif][old_motif] > 40:
                if new_motif.split(" ")[0] in allowed_motifs and old_motif.split(" ")[0] in allowed_motifs:
                    unique_motifs.add(old_motif.split("_")[0])
                    unique_motifs.add(new_motif.split("_")[0])
                    G.add_edge(
                        new_motif,
                        old_motif,
                        weight=graphing_obj[new_motif][old_motif]
                        / max(motif_appears_in_replicates[new_motif], motif_appears_in_replicates[old_motif]),
                    )
    print(unique_motifs)
    # Extract weights for normalization
    weights = [d["weight"] for (u, v, d) in G.edges(data=True)]
    weights = np.sqrt(weights)
    min_alpha, max_alpha = 0.2, 1.0
    norm_weights = [
        (weight - min(weights)) / (max(weights) - min(weights)) * (max_alpha - min_alpha) + min_alpha
        for weight in weights
    ]

    # Extract node scores for normalization
    scores = [node_scores[node] for node in G.nodes()]
    norm_scores = (np.array(scores) - min(scores)) / (max(scores) - min(scores))

    # Generate colors based on normalized scores
    cmap = cm.get_cmap("Spectral_r")
    colors = cmap(norm_scores)

    # Determine node sizes based on motif sizes
    sizes = [motif_appears_in_replicates[node] for node in G.nodes()]
    min_size, max_size = 100, 1000  # You can adjust these values as needed
    norm_sizes = [(size - min(sizes)) / (max(sizes) - min(sizes)) * (max_size - min_size) + min_size for size in sizes]

    # Map nodes to their colors for edge coloring
    node_color_map = {node: colors[i] for i, node in enumerate(G.nodes())}

    # Draw the graph
    pos = nx.spring_layout(G, seed=7)  # positions for all nodes - seed for reproducibility

    # Adjust positions based on generations
    gen_values = [generations[node] for node in G.nodes()]
    print(min(gen_values), max(gen_values))
    norm_gens = (np.array(gen_values) - min(gen_values)) / (max(gen_values) - min(gen_values))
    print(len(pos))
    for node in pos:
        # print(node)
        ord_num = ord(node[0]) - ord("A")
        print(ord_num, node[0])
        if ord_num == 2:
            ord_num = random.randint(-5, 0)
        if ord_num == 11:  # L
            ord_num = random.randint(1, 5)
        if ord_num == 13:  # N
            ord_num = random.randint(7, 11)
        if ord_num == 15:  # P
            ord_num = random.randint(13, 30)
        if ord_num == 18:  # S
            ord_num = random.randint(37, 44)
        if ord_num == 25:
            ord_num = random.randint(48, 52)
        """
        if ord_num > 8 and ord_num < 18:
            ord_num += random.randint(-20, 20) / 3
        else:
            ord_num += random.randint(-10, 10)
        """
        pos[node][1] = ord_num
        pos[node][0] = norm_gens[list(G.nodes()).index(node)]

    sorted_nodes = sorted(G.nodes(), key=lambda node: pos[node][0])
    # for sorted_nodes
    nx.draw_networkx_nodes(G, pos, node_size=norm_sizes, node_color=colors, cmap="viridis")
    for (u, v, d), alpha in zip(G.edges(data=True), norm_weights):
        nx.draw_networkx_edges(
            G, pos, edgelist=[(u, v)], width=int(alpha * 6), alpha=alpha, edge_color=[node_color_map[u]]
        )

    nx.draw_networkx_labels(G, pos, font_size=5, font_family="sans-serif")

    # Draw edge labels
    edge_labels = nx.get_edge_attributes(G, "weight")
    # nx.draw_networkx_edge_labels(G, pos, edge_labels, font_size=4)
    sm = plt.cm.ScalarMappable(cmap="Spectral_r", norm=plt.Normalize(vmin=min(scores), vmax=max(scores)))
    sm.set_array([])  # Only needed for ScalarMappable
    plt.colorbar(sm, ax=plt.gca(), label="Node Scores")
    # Display the plot
    ax = plt.gca()
    ax.margins(0.08)
    plt.axis("off")
    plt.tight_layout()
    plt.show()


bed_lines = []
gen_add = 0
survivors = 4
line_seq_holder = {}


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
with open("fimo_40_160/fimo.tsv") as file:
    for line in file:
        cell = line.split("\t")
        if cell[0] != "motif_id" and len(cell) > 5:
            i_line, gen, surv_rep, score = cell[2].split("_")
            # chr_n = ord(i_line[-1]) - ord("a") + 1
            chr_n = i_line.split("rep")[-1]
            chrom = f"chr{chr_n}"
            gen_n = int(gen.split("G")[-1])
            # print(cell)
            if float(cell[-3]) < 10**-6:
                # print(cell[0].split("_")[0])
                cell[0] = cell[0].split("_")[0] + "_" + score
                motif, score = cell[0].split("_")
                b_line = "\t".join([chrom, cell[3], cell[4], cell[0], "1", cell[5]]) + "\n"
                line_name = i_line + "_" + surv_rep
                m_pos = int((int(cell[3]) + int(cell[4])) / 2)
                line_seq_holder[line_name][gen].append((motif, m_pos))
                # print(line_seq_holder[line_name][gen])
                # bed_lines.append(b_line)

violin_dict = {}
color_motif = {}
clrs = sns.color_palette("deep", 10) + sns.color_palette("husl", 191)
clr_i = 0
graphing_obj = {}
node_scores = {}
generations = {}
motif_appears_in_replicates = {}
for line in line_seq_holder:
    old_motifs = set()
    for gen in line_seq_holder[line]:
        gen_num = int(gen.replace("G", ""))
        last_gen = f"G{gen_num-1}"
        if gen != "G1":
            try:
                print(gen)
                score = math.sqrt(float(line_seq_holder[line][gen][0][1]))
                last_score = math.sqrt(float(line_seq_holder[line][last_gen][0][1]))
                seq = line_seq_holder[line][gen][0][0]
                motifs = line_seq_holder[line][gen][1:]

                # print(gen, score, motifs)
                motif_holder = {}
                if len(motifs) == 0:
                    motifs.append(("None", 0))
                print(gen, len(motifs), motifs)
                for motif, m_pos in motifs:
                    if motif not in motif_holder:
                        motif_holder[motif] = 0
                    motif_holder[motif] += 1
                motifs = set()
                for motif in motif_holder:
                    motif_id = f"{motif} ({motif_holder[motif]})"
                    if motif_id not in motif_appears_in_replicates:
                        motif_appears_in_replicates[motif_id] = set()
                    motif_appears_in_replicates[motif_id].add(line)

                    """
                    if motif_id not in node_scores:
                        node_scores[motif_id] = []
                    node_scores[motif_id].append(score)
                    """
                    motifs.add(motif_id)
                if motifs != old_motifs:
                    new_motifs = motifs - old_motifs
                    if len(new_motifs) != 0:
                        # for motif in new_motifs:
                        for motif in new_motifs:
                            print(motif)
                            if motif not in node_scores:
                                node_scores[motif] = []
                            node_scores[motif].append(score)
                            if motif not in generations:
                                generations[motif] = []
                            generations[motif].append(gen_num)
                            if motif not in graphing_obj:
                                graphing_obj[motif] = {}
                            for o_motif in old_motifs:
                                if o_motif not in graphing_obj[motif]:
                                    graphing_obj[motif][o_motif] = 0
                                graphing_obj[motif][o_motif] += 1
                        # print(gen, old_motifs, motifs)

                old_motifs = motifs
            except Exception:
                pass
print(graphing_obj)
motif_sizes = {}
for node in node_scores:
    if node not in generations:
        generations[node] = -1
    generations[node] = np.median(generations[node])
    node_scores[node] = np.median(node_scores[node])
    motif_sizes[node] = (int(node.split("(")[-1].split(")")[0]) * 100) ** 3
    motif_appears_in_replicates[node] = len(motif_appears_in_replicates[node])

net_x(graphing_obj, node_scores, motif_sizes, generations, motif_appears_in_replicates)
