import os
import numpy as np
import pyBigWig
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob
import re

# Define main regions of interest
regions = [("chr1", 2448, 3481), ("chr1", 2448, 2964), ("chr1", 2964, 3482)]
# regions = [("chr1", 2448, 2975), ("chr1", 2448, 2448 + 218), ("chr1", 2448 + 218 + 89, 2975)]  # zhh
region_labels = ["Both Enhancers", "Enhancer 1", "Enhancer 2"]

# Define control regions
control_regions = [("chr1", 1112, 1343), ("chr1", 3890, 4120)]
# control_regions = [("chr1", 1112, 1343), ("chr1", 3384, 3614)]  # shh

# Directory and file selection
bigwig_dir = "../data/rep6 del_supp"
# bigwig_dir = "/Users/blake/Documents/gargLab/onescreen/bws_onescreen_atacrep6wdel"
# if ("g90e" in f.lower()) and ("r6g90align" not in f.lower()) and "313" not in f.lower() or "r6g1align" in f.lower()
# bigwig_dir = "/Users/blake/Documents/gargLab/onescreen/bws_onescreen_atacrep6wdel_2"
# bw_files = [f for f in bw_files_messy if ("" in f.lower()) and ("r6g90" not in f.lower()) and "313" not in f.lower()]

# r6g1align_name_1, r28g1
bw_files_messy = glob(os.path.join(bigwig_dir, "*.bw"))

# Filter relevant files
bw_files = [f for f in bw_files_messy if ("" in f.lower())]


def extract_g_number(filename):
    match = re.search(r"_G(\d+)", filename)
    if match:
        pass
        # if int(match.group(1)) == 90:
        #    return -1
    return int(match.group(1)) if match else float("inf")


# Sort files by G number
bw_files = sorted(bw_files)
bw_files = sorted(bw_files, key=extract_g_number)

bw_files = [
    "../data/rep6 del_supp/r6g90align_cpm_name_1.SS_6_1_G90.bw",
    "../data/rep6 del_supp/osn1align_cpm_name_1.SS_OSN1Del.bw",
    "../data/rep6 del_supp/osn2delalign_cpm_name_1.SS_OSN2Del.bw",
    "../data/rep6 del_supp/osn3delalign_cpm_name_1.SS_OSN3Del.bw",
    "../data/rep6 del_supp/zic3delalign_cpm_name_1.SS_6_90_ZIC3Del.bw",
    "../data/rep6 del_supp/fourdelalign_cpm_name_1.SS_4Del.bw",
    "../data/rep6 del_supp/gcdelalign_cpm_name_1.SS_GCDel.bw",
    "../data/rep6 del_supp/fourgcdelalign_cpm_name_1.SS_4GCDel.bw",
    "../data/rep6 del_supp/r6g1align_cpm_name_1.SS_6_1_G1.bw",
]
# Set reference file
ref_file = next(f for f in bw_files_messy if "r6g90align" in f.lower())


def get_region_values(bw_path, chrom, start, end):
    with pyBigWig.open(bw_path) as bw:
        values = bw.values(chrom, start, end)
        return np.array([0 if np.isnan(v) else v for v in values])


# --- NEW: Get mean control value per sample ---
def get_mean_control_value(bw_path):
    values = []
    for region in control_regions:
        region_vals = get_region_values(bw_path, *region)
        values.append(np.mean(region_vals))
    return np.mean(values)


# Compute values for all regions and normalize by control
data_by_region = []
ref_data_by_region = []

# Precompute control means
control_means = [get_mean_control_value(bw) for bw in bw_files]
ref_control_mean = get_mean_control_value(ref_file)
print(control_means)

for region in regions:
    data = [get_region_values(bw, *region) for bw in bw_files]
    ref_data = get_region_values(ref_file, *region)
    data_by_region.append(data)
    ref_data_by_region.append(ref_data)

labels = [os.path.basename(f).split("_")[-1].replace(".bw", "") for f in bw_files]


# Prepare matrix: rows = samples, columns = regions
heatmap_matrix = []

for i in range(len(bw_files)):
    row = []
    for j in range(len(regions)):
        sample_data = data_by_region[j][i]
        ref_data = ref_data_by_region[j]
        epsilon = 1e-1
        # norm_sample_data = sample_data / (control_means[i] + epsilon)
        # norm_ref_data = ref_data / (ref_control_mean + epsilon)
        ratio = (sample_data + epsilon) / (ref_data + epsilon)
        log2_diff = np.log2(ratio)
        row.append(np.mean(log2_diff))
    heatmap_matrix.append(row)

# Convert to numpy array
heatmap_matrix = np.array(heatmap_matrix)

# Plot single heatmap
plt.figure(figsize=(4.8, len(bw_files) * 0.5 + 2))
sns.heatmap(
    heatmap_matrix,
    yticklabels=labels,
    xticklabels=region_labels,
    center=0.0,
    vmax=3,
    vmin=-3,
    cmap="coolwarm",
    annot=True,
    fmt=".2f",
    cbar=True,
)
# plt.title(f"log2 fold-change vs {os.path.basename(ref_file)}")
# plt.ylabel("Samples")
# plt.xlabel("Region")
plt.tight_layout()
plt.savefig("../figures/suppdel.svg")
plt.show()
