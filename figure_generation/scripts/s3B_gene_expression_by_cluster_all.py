from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats

ce_to_cluster = {}

cluster_labels = {}


def get_gene_connection():
    gene_expression_dict = get_expression_count()
    # print(gene_expression_dict)
    ce_to_expression = {}
    with open("../data/ce_gene_closest.bed") as file:
        for line in file:
            # print(line)
            cell = line.strip().split("\t")
            ce = f"{cell[0]}:{cell[1]}-{cell[2]}"
            # print(cell)
            if ce not in ce_to_expression:
                ce_to_expression[ce] = []
            for gene in [cell[-2]]:
                gene = gene.lower()
                if gene in gene_expression_dict:
                    ce_to_expression[ce].append((gene, gene_expression_dict[gene]))
    return ce_to_expression
    # print(ce_to_expression)


def get_expression_count():
    gene_expression_dict = {}
    # data/whyte2013_expression.txt
    # "../FIRE_story_v1/data/esc_expression_level_3states.tsv"
    with open("../data/esc_expression_level_3states.tsv") as file:
        cell_num = -1
        for line in file:
            print(line)
            cell = line.split("\t")
            cell[cell_num] = cell[cell_num].strip()

            if cell[cell_num] == "" or cell[cell_num] == "NA":
                cell[cell_num] = float(0)
            else:
                cell[cell_num] = float(cell[cell_num])
            """
            if cell[6] > 10:
                print(cell[0].lower())
            """
            if cell[0].lower() not in gene_expression_dict:
                gene_expression_dict[cell[0].lower()] = cell[cell_num]
            gene_expression_dict[cell[0].lower()] = max(gene_expression_dict[cell[0].lower()], cell[cell_num])

    return gene_expression_dict


def get_CEs_by_cluster_and_gene_expression():
    data_holder = {}
    cluster_labels = {}
    fig = plt.figure(figsize=(7, 3))
    ax = fig.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot
    ax.spines[["right", "top"]].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.5)  # Bottom (x-axis) spine
    ax.spines["left"].set_linewidth(1.5)  # Left (y-axis) spine
    ax.tick_params(axis="both", width=2)
    ce_to_expression = get_gene_connection()
    # print(ce_to_expression)
    with open("../data/fire_3runs_BOX.svg3.bed") as file:
        for line in file:
            if "deepTools_group" not in line:
                cell = line.strip().split("\t")
                ce_name = f"{cell[0]}:{cell[1]}-{cell[2]}"
                cluster_labels[ce_name] = int(cell[-1][-1])
    used_genes = set()
    used_stitches = set()
    with open("../data/3runs_ce_rank.txt") as file:
        for line in file:
            if ":" in line:
                ceA, ceB, stitch, OR, rank, isSuper = line.split("\t")
                OR = float(OR)
                if isSuper.strip() != "Super":
                    isSuper = "w"
                c_id = "_".join(map(str, sorted([cluster_labels[ceA], cluster_labels[ceB]]))) + f"{isSuper.strip()}"
                c_id_no_super = "_".join(map(str, sorted([cluster_labels[ceA], cluster_labels[ceB]])))
                cell = line.strip().split("\t")
                # print(c_id, OR, isSuper, ce_to_expression[ceA], ce_to_expression[ceB])
                if c_id not in data_holder:
                    data_holder[c_id] = []
                # print(ceA, ceB)
                # print(ce_to_expression)
                maxA = max(ce_to_expression[ceA] + [("none", 0)], key=lambda x: x[1])
                maxB = max(ce_to_expression[ceB] + [("none", 0)], key=lambda x: x[1])
                if stitch + c_id_no_super not in used_stitches:
                    used_stitches.add(stitch + c_id_no_super)
                    if maxA[0] != "none":
                        if maxA[0] + c_id_no_super not in used_genes:
                            data_holder[c_id].append(maxA[1])
                            # used_genes.add(maxA[0] + c_id_no_super)
                    else:
                        pass  # data_holder[c_id].append(0)
                    if maxB[0] != "none":
                        if maxB[0] + c_id_no_super not in used_genes:
                            data_holder[c_id].append(maxB[1])
                            # used_genes.add(maxB[0] + c_id_no_super)
                    else:
                        pass  # data_holder[c_id].append(0)
    # print(used_genes)
    # data_holder[c_id].append(ce_to_expression[ceA])
    # data_holder[c_id].append(ce_to_expression[ceB])
    box_clusters = []
    box_labels = []
    data_holder = dict(sorted(data_holder.items()))
    print(data_holder["1_1Super"])
    # print(data_holder["1_1"])
    # print(stats.ttest_ind(data_holder["1_1Super"], data_holder["1_1"]))
    for name in data_holder:
        # if "1_1" in name or "2_2" in name or "3_3" in name:
        if "_" in name:
            print(name, len(data_holder[name]))
            if "Super" not in name:
                arrs = np.array_split(data_holder[name], 1)
                for i, a in enumerate(arrs):
                    box_clusters.append(a)
                    name = name.replace("_", ",").replace("w", " ORE")
                    box_labels.append(f"{name}")
            else:
                box_clusters.append(data_holder[name])
                name = name.replace("_", ",").replace("Super", " HCRE")
                box_labels.append(name)
    # plt.yscale("log")
    """
    # rgb(206,70,40)
    --black: #010000ff;
    --dark-cyan: #189088ff;
    --persian-red: #D0373Aff;
    --earth-yellow: #E4AA56ff;
    --silver: #BCBCBDff;
    """
    sns.set(font_scale=4)
    sns.set_theme(style="whitegrid")
    # custom_colors = ["#de3e2b", "#BCBCBD"]
    colors = ["#777777", "#d9d9d9"] * 20
    ax = sns.boxplot(
        data=box_clusters,
        width=0.5,
        linewidth=1.2,
        fliersize=0,
        widths=0.25,
        palette=colors,
        boxprops={"edgecolor": "black"},  # Box outline color
        whiskerprops={"color": "black"},  # Whisker color
        capprops={"color": "black"},  # Cap color
        medianprops={"color": "black"},  # Median line color
    )
    plt.ylim(([-999, 7000]))
    ax.set_xticklabels(box_labels, rotation=-30, ha="center")
    plt.tight_layout()
    plt.savefig("../figures/s3b_cluster_all.svg")
    plt.show()


# get_gene_connection()
get_CEs_by_cluster_and_gene_expression()
