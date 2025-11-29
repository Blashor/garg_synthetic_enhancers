import seaborn as sns
from matplotlib import pyplot as plt
import math
import numpy

enh_len = []
enh_spacing = []
super_stuff = {}
with open("../FIRE_story_v1/data/ce_rank.txt") as file:
    for line in file:
        if "Super" in line:
            cells = line.split("\t")
            ce_coords = []
            for ce in cells[:2]:
                start, end = map(int, ce.split(":")[1].split("-"))
                ce_coords.append((start, end))
                enh_len.append(end - start)
            spacing = ce_coords[1][0] - ce_coords[0][1]
            if spacing > 0:
                print(cells[2], cells[3], spacing)
                enh_spacing.append(spacing)


print(numpy.median(enh_len))
print(numpy.median(enh_spacing))
sns.ecdfplot(data=enh_len)
plt.show()

sns.ecdfplot(data=enh_spacing)
plt.show()
