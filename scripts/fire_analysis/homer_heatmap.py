import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import linkage
from matplotlib.colors import ListedColormap

# File path
file_path = "/Users/blake/Downloads/3runs-FDR-FIRE-peaks.bed.txt"
# file_path = "/Users/blake/Downloads/annotated_super_ce_motifs_500.txt"
# Read file and extract motif names
with open(file_path) as file:
    lines = file.readlines()

# Extract motif names from the first row
motif_names = list(map(lambda x: x.split("(")[0], lines[0].strip().split("\t")[21:]))

# Identify indices of motifs that match the filter list
filter_list = []
with open("/Users/blake/Downloads/knownResults.txt") as file:
    for i, line in enumerate(file):
        cell = line.split("(")[0]
        filter_list.append(cell)
        if i == 50:
            break
filter_motifs = set(filter_list)

filtered_indices = [i for i, name in enumerate(motif_names) if name in filter_motifs]
filtered_motif_names = [motif_names[i] for i in filtered_indices]

# Process the data into a NumPy array
fire_homer = []
for line in lines[1:]:  # Skip header
    cells = line.split("\t")
    fire_homer.append([x.count("(") if x.strip() != "" else 0 for x in cells[21:]])


# Convert to NumPy array and filter by selected motifs
fire_homer_array = np.array(fire_homer, dtype=int)[:, filtered_indices]

# Create a pandas DataFrame
df = pd.DataFrame(fire_homer_array, columns=filtered_motif_names)

# KMeans clustering for rows
kmeans = KMeans(n_clusters=7, random_state=42)  # Adjust n_clusters as needed
row_clusters = kmeans.fit_predict(df)

# Create a color map for the KMeans clusters
cmap = ListedColormap(sns.color_palette("Set2", n_colors=7).as_hex())  # Change 'Set2' to a different palette if needed
row_colors = [cmap(i) for i in row_clusters]

# Reorder rows based on KMeans clustering (so similar rows are together)
df_sorted = df.iloc[row_clusters.argsort()]  # Sort rows based on KMeans clustering

# Reorder row_colors based on the same sorting
row_colors_sorted = [row_colors[i] for i in row_clusters.argsort()]

# Hierarchical clustering for columns (linkage)
column_linkage = linkage(df.T, method="ward", metric="euclidean", optimal_ordering=True)

# Create a heatmap using clustermap only for columns (row clustering is handled by KMeans)
plt.figure(figsize=(12, 8))
sns.clustermap(
    df_sorted,
    row_cluster=False,  # Disable row clustering (KMeans already handled)
    col_cluster=True,  # Enable column clustering
    row_colors=row_colors_sorted,  # Color rows based on KMeans clustering (sorted)
    figsize=(12, 8),
    cmap="Reds",
    vmax=5,
    xticklabels=True,
    yticklabels=False,  # Disable yticklabels for better readability
    col_linkage=column_linkage,  # Use column linkage for hierarchical clustering
)

# Save the heatmap
plt.title("Heatmap with KMeans for Rows and Column Linkage Clustering")
plt.xlabel("Motifs")
plt.ylabel("Samples")
plt.savefig("/Users/blake/Downloads/fire_homer_clustermap_columns_sorted.png", dpi=300)
plt.show()
