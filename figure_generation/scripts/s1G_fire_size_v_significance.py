import matplotlib as mpl
import seaborn as sns
import math
import numpy
import mpl_style


mpl.use("Agg")
from matplotlib import pyplot as plt


def main(input_bed="../data/3runs-FIRE.bed", output="../figures/s2e_3runs-FIRE.svg"):
    print("go")
    # plt.figure()
    mpl_style.set_style()
    x = []
    y = []
    window_size = 25
    bin_obj = {}
    with open(input_bed) as file:
        for line in file:
            cell = line.split("\t")
            if cell[-2] != "1.0":
                try:
                    FIRE_score = float(cell[-2])
                    if FIRE_score <= 0.1:
                        length = int(cell[2]) - int(cell[1])
                        if length < 1500:
                            bin_itr = math.floor(length / window_size)
                            if bin_itr not in bin_obj:
                                bin_obj[bin_itr] = []
                            bin_obj[bin_itr].append(FIRE_score)
                except ValueError:
                    pass

    print("Done")
    myKeys = list(bin_obj.keys())
    print(len(myKeys))
    myKeys.sort()
    bin_obj = {i: bin_obj[i] for i in myKeys}

    for bin_itr in bin_obj:
        x.append(window_size * bin_itr)
        print(
            window_size * bin_itr,
            len(bin_obj[bin_itr]),
        )
        y.append(bin_obj[bin_itr][:10000])

    ax = sns.boxplot(data=y, fliersize=0, linewidth=0.2, color=(157 / 255, 30 / 255, 30 / 255), width=0.8)

    ax.set_xticklabels(x)
    for ind, label in enumerate(ax.get_xticklabels()):
        if ind % 10 == 0:  # every 10th label is kept
            label.set_visible(True)
        else:
            label.set_visible(False)
    plt.xlabel("Element Length")
    plt.ylabel("FIRE Confidence Score")
    plt.savefig(output, dpi=300)


main()
