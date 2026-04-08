from matplotlib import pyplot as plt
import seaborn as sns
from scipy.stats.stats import pearsonr

from scipy.optimize import curve_fit
import math
from matplotlib.colors import LinearSegmentedColormap

ys = []
xs = []
with open("../data/t1_ce_43_scatter.txt") as file:
    for line in file:
        cell = line.split(" ")
        # print(cell)
        if len(cell) == 2:
            # print(cell[0])
            x = float(cell[0])
            y = float(cell[1])
            xs.append(x)
            ys.append(y)
# sns.scatterplot(x=xs, y=ys, linewidth=0, s=5, alpha=0.5)
# plt.show()

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# xs = 2 ** np.array(xs)
# ys = 2 ** np.array(ys)
xs = np.array(xs)
ys = np.array(ys)
# Calculate line of best fit and Pearson correlation
slope, intercept, r_value, p_value, std_err = linregress(xs, ys)
from scipy.stats import spearmanr

# xs and ys are your data arrays/lists
corr, p_value = pearsonr(xs, ys)
print(r_value, corr, p_value)
# Create the line of best fit
line = slope * xs + intercept

# Plot the points
sample = 10000
sample_small = 1500
indexes = np.random.choice(np.arange(0, len(xs)), sample, replace=False)
indexes2 = np.random.choice(np.arange(0, len(xs)), sample_small, replace=False)
# print(indexes)
xs_small = xs[indexes]

xs_v_small = xs[indexes2]
ys_v_small = ys[indexes2]
ys = ys[indexes]
fig = plt.figure(figsize=(3.5, 3.5))

ax = fig.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot
ax.spines[["right", "top"]].set_visible(False)
ax.spines["bottom"].set_linewidth(1.5)  # Bottom (x-axis) spine
ax.spines["left"].set_linewidth(1.5)  # Left (y-axis) spine
ax.tick_params(axis="both", width=2)
gold = sns.color_palette("husl", 8)[1]
indexes = xs_small != 0
indexes2 = xs_v_small != 0

xs_v_small = xs_v_small[indexes2]
ys_v_small = ys_v_small[indexes2]
xs_small_filter = xs_small[indexes]
ys_small = ys[indexes]

sns.kdeplot(
    x=xs_small_filter,
    y=ys_small,
    cmap=LinearSegmentedColormap.from_list("custom_cmap", ["white", (157 / 255, 36 / 255, 36 / 255)]),
)

sns.scatterplot(
    x=xs_v_small,
    y=ys_v_small,
    s=2,
    alpha=0.25,
    color=(157 / 255, 36 / 255, 36 / 255),
    edgecolor=(157 / 255, 36 / 255, 36 / 255),
)

ax = sns.regplot(
    x=xs_small,
    y=ys,
    color=None,
    order=1,
    scatter_kws={"color": None, "s": 10, "alpha": 0},
    line_kws={"linestyle": "--", "color": (157 / 255, 36 / 255, 36 / 255)},
)


# plt.scatter(xs_small, ys, color="grey", label="Data Points", s=1, alpha=0.5)

# Plot the line of best fit
# plt.plot(xs, line, color="red", label=f"Line of Best Fit\n$R^2={r_value**2:.2f}$")

# Add labels, legend, and title
plt.xlabel("Log2 Prediction")
plt.ylabel("Log2 Experiment")
# plt.title("Line of Best Fit and Data Points")
# plt.legend()
# plt.ylim([0, 8])
# plt.xlim([0, 8])

# Show the plot
plt.savefig("../figures/2c_scatter.svg")
plt.show()

# Print Pearson correlation coefficient
print(f"Pearson correlation coefficient (R): {r_value}")
print(f"R-squared value: {r_value**2}")
