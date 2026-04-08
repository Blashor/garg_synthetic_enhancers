import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import os
import seaborn as sns
import math
import matplotlib.cm as cm
import random


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

allowed_motifs = set(["SOX2", "POU5F1", "NANOG", "BACH2", "CTCF", "ZIC3", "LEF1", "POU2F2", "POU3F1", "POU2F1", "None"])


def net_x(graphing_obj, node_scores, motif_sizes, generations, motif_appears_in_replicates):
    unique_motifs = set()
    G = nx.Graph()

    # Add edges based on the graphing object
    for new_motif in graphing_obj:
        for old_motif in graphing_obj[new_motif]:
            if node_scores[new_motif] < 1500000 and node_scores[old_motif] < 15000000:
                # graphing_obj[new_motif][old_motif] > 40 or
                # if graphing_obj[new_motif][old_motif] > 40:
                if new_motif.split(" ")[0] in allowed_motifs and old_motif.split(" ")[0] in allowed_motifs:
                    unique_motifs.add(old_motif.split("_")[0])
                    unique_motifs.add(new_motif.split("_")[0])
                    G.add_edge(
                        new_motif,
                        old_motif,
                        weight=np.mean([graphing_obj[new_motif][old_motif]]),
                    )
    print(unique_motifs)
    # Extract weights for normalization
    weights = [d["weight"] for (u, v, d) in G.edges(data=True)]
    # weights = np.sqrt(weights)

    min_alpha, max_alpha = 0.1, 1.0
    # sys.exit()
    norm_weights = [(weight) / (160) * (max_alpha - min_alpha) + min_alpha for weight in weights]

    # Extract node scores for normalization
    scores = [node_scores[node] for node in G.nodes()]
    norm_scores = (np.array(scores) - min(scores)) / (max(scores) - min(scores))

    # Generate colors based on normalized scores
    motif_cmap = sns.color_palette() + [(0, 0, 0, 0.9)]
    cmap = cm.get_cmap("Spectral_r")
    colors = cmap(norm_scores)

    list_unique = ["SOX2", "POU5F1", "NANOG", "BACH2", "CTCF", "ZIC3", "LEF1", "POU2F2", "POU3F1", "POU2F1", "None"]
    palette_big = sns.color_palette("hls", 24)
    palette_small = sns.color_palette("hls", 7)
    motif_color_map = {
        "POU5F1": palette_big[0],
        "SOX2": palette_big[1],
        "NANOG": palette_big[2],
        "BACH2": palette_big[18],
        "LEF1": palette_big[16],
        "ZIC3": palette_small[2],
        "CTCF": palette_small[4],
        "POU2F2": palette_big[11],
        "POU3F1": palette_big[12],
        "POU2F1": palette_big[13],
        "None": (0, 0, 0, 0.9),
    }

    motif_node_colors = []
    # print(colors)
    # sys.exit()
    for node in G.nodes():
        motif_node_colors.append(motif_color_map[node.split(" ")[0]])
    # Determine node sizes based on motif sizes
    sizes = [motif_appears_in_replicates[node] for node in G.nodes()]
    min_size, max_size = 100, 1000  # You can adjust these values as needed
    norm_sizes = [(size - min(sizes)) / (max(sizes) - min(sizes)) * (max_size - min_size) + min_size for size in sizes]
    print(norm_sizes)
    # Map nodes to their colors for edge coloring
    node_color_map = {node: motif_color_map[node.split(" ")[0]] for i, node in enumerate(G.nodes())}
    # node_color_map = {node: colors[i] for i, node in enumerate(G.nodes())}
    # Draw the graph
    pos = nx.spring_layout(G, seed=7)  # positions for all nodes - seed for reproducibility

    # Adjust positions based on generations
    gen_values = [generations[node] for node in G.nodes()]
    print(min(gen_values), max(gen_values))
    norm_gens = (np.array(gen_values) - min([0])) / (max(gen_values) - min([0]))
    print(len(pos))
    node_pos = {
        "SOX2": 1 - 0.5,
        "POU5F1": 2,
        "NANOG": 3 + 0.5,
        "CTCF": -3,
        "ZIC3": -2,
        "LEF1": -1,
        "POU2F2": 5,
        "POU3F1": 6,
        "POU2F1": 7,
        "BACH2": 8,
        "None": 5,
    }
    for node in pos:
        # print(node)
        node_select = node.split(" ")[0]
        ord_num = node_pos[node_select] * (1 + (random.randint(0, 30) / 100))  # ord(node[0]) - ord("A")
        if node == "CTCF (1)":
            ord_num = -2.5
        if node == "CTCF (2)":
            ord_num = -3.5
        if node == "CTCF (5)":
            ord_num = -3
        if node == "CTCF (3)":
            ord_num = -3.7
        if node == "LEF1 (1)":
            ord_num = -0.6
        if node == "CTCF (4)":
            ord_num = -4.6
        if node == "POU5F1 (6)":
            ord_num = 1.5
        if node == "POU5F1 (5)":
            ord_num = 2.5
        if node == "ZIC3 (1)":
            ord_num = -1.9
        if node == "ZIC3 (3)":
            ord_num = -2.8
        if node == "ZIC3 (4)":
            ord_num = -2.0
        if node == "SOX2 (5)":
            ord_num = 1
        if node == "NANOG (1)":
            ord_num = 5
        """
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
    print(motif_node_colors)
    nx.draw_networkx_nodes(
        G, pos, node_size=norm_sizes, node_color=colors, edgecolors=motif_node_colors, linewidths=3, cmap="viridis"
    )
    for (u, v, d), alpha in zip(G.edges(data=True), norm_weights):
        mixed_color = (np.array(node_color_map[u])[:3] + np.array(node_color_map[v])[:3]) / 2
        if node_color_map[u][0] == 0:
            mixed_color = (0, 0, 0)
        print(alpha)
        nx.draw_networkx_edges(
            G, pos, edgelist=[(u, v)], width=int(alpha * 10), alpha=min([0.8, alpha]), edge_color=[mixed_color]
        )

    labels = {node: node.split(" ")[1].strip("()") for node in G.nodes()}

    nx.draw_networkx_labels(G, pos, labels=labels, font_size=13, font_family="sans-serif")

    # Draw edge labels
    edge_labels = nx.get_edge_attributes(G, "weight")
    # nx.draw_networkx_edge_labels(G, pos, edge_labels, font_size=4)
    sm = plt.cm.ScalarMappable(cmap="Spectral_r", norm=plt.Normalize(vmin=min(scores), vmax=max(scores)))
    sm.set_array([])  # Only needed for ScalarMappable
    plt.colorbar(sm, ax=plt.gca(), label="Odds Ratio")
    # Display the plot
    ax = plt.gca()
    ax.spines[["right", "top"]].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.5)
    ax.spines["left"].set_linewidth(1.5)
    ax.tick_params(axis="both", width=2)
    ax.margins(0.08)
    plt.xticks(np.arange(0, 1, 0.2))
    # plt.axis("off")
    plt.tight_layout()
    plt.show()


bed_lines = []
gen_add = 0
survivors = 4
line_seq_holder = {}


tsv_lines = []
# INSIDE_256_g117.out
with open("../data/INSIDE_jun13.txt") as file:
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
        if int(cell[0].split("Gen:")[-1]) < 160:
            cell[0] = "G" + str(int(cell[0].split("Gen:")[-1]) + gen_add)
            current_stat += float(cell[3])
            cell[1] += "_" + str(surv_num)
            if cell[1] not in line_seq_holder:
                line_seq_holder[cell[1]] = {}
            if cell[0] not in line_seq_holder[cell[1]]:
                line_seq_holder[cell[1]][cell[0]] = []
            line_seq_holder[cell[1]][cell[0]].append((cell[2], cell[3]))
# print(line_seq_holder)
with open("../data/INSIDE_40_definitive/fimo.tsv") as file:
    for line in file:
        cell = line.split("\t")
        if cell[0] != "motif_id" and len(cell) > 5:
            i_line, gen, surv_rep, score = cell[2].split("_")
            # chr_n = ord(i_line[-1]) - ord("a") + 1
            chr_n = i_line.split("rep")[-1]
            chrom = f"chr{chr_n}"
            gen_n = int(gen.split("G")[-1])
            # print(cell)
            if float(cell[-3]) < 10**-6 and gen_n < 160:
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
                    motif = motif.replace("PO", "POU")
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
    print(node, len(motif_appears_in_replicates[node]))
    motif_appears_in_replicates[node] = len(motif_appears_in_replicates[node])

net_x(graphing_obj, node_scores, motif_sizes, generations, motif_appears_in_replicates)
