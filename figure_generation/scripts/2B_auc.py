from matplotlib import pyplot as plt
import seaborn as sns
import math
import numpy as np

fpr = np.loadtxt("../data/1_fpr.txt", delimiter=",", skiprows=1)  # Skip the header row
tpr = np.loadtxt("../data/1_tpr.txt", delimiter=",", skiprows=1)

fig = plt.figure(figsize=(3, 3))
plt.tight_layout()
ax = fig.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot
ax.spines[["right", "top"]].set_visible(False)
ax.spines["bottom"].set_linewidth(1.5)  # Bottom (x-axis) spine
ax.spines["left"].set_linewidth(1.5)  # Left (y-axis) spine
ax.tick_params(axis="both", width=2)

plt.plot([0, 1], [0, 1], c="black", linewidth=1.5, linestyle="--", alpha=0.7)
plt.plot(fpr, tpr, c="black")
ax = plt.gca()
ax.set_xlabel("False positive rate")
ax.set_ylabel("True positive rate")
ax.text(0.99, 0.02, "AUROC %.3f" % 0.901, horizontalalignment="right")  # , fontsize=14)
plt.savefig("../figures/2b_auc.svg")
plt.show()
