############################################################
# BREAST CANCER scRNA INTERPRETATION ANALYSIS (UPDATED)
############################################################

import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

############################################################
# OUTPUT DIRECTORIES
############################################################

BASE = "scRNA_results/interpretation"

dirs = [
f"{BASE}/plots",
f"{BASE}/tables"
]

for d in dirs:
    os.makedirs(d, exist_ok=True)

print("Output folders created")

############################################################
# LOAD DATASET
############################################################

adata = sc.read_h5ad("scRNA_results/data/BC_data_celltypist.h5ad")

print("Dataset loaded:")
print(adata)

############################################################
# CELL TYPE COMPOSITION
############################################################

print("\nRunning cell composition analysis")

cell_counts = pd.crosstab(
adata.obs["majority_voting"],
adata.obs["condition"]
)

cell_props = cell_counts.div(
cell_counts.sum(axis=0),
axis=1
) * 100

cell_props.to_excel(
f"{BASE}/tables/celltype_proportions.xlsx"
)

print("Cell proportions saved")

############################################################
# CELL TYPE COMPOSITION PLOT
############################################################

plt.figure(figsize=(12,7))

cell_props.T.plot(
kind="bar",
stacked=True
)

plt.ylabel("Cell proportion (%)")
plt.title("Cell Type Composition Across Conditions")

plt.legend(
bbox_to_anchor=(1.02,1),
loc="upper left"
)

plt.tight_layout()

plt.savefig(
f"{BASE}/plots/celltype_composition.png"
)

plt.close()

############################################################
# BRCA1 vs BRCA2 DIFFERENTIAL EXPRESSION
############################################################

print("\nRunning BRCA1 vs BRCA2 DEG analysis")

tumor = adata[
adata.obs["condition"].isin(["BRCA1","BRCA2"])
].copy()

sc.tl.rank_genes_groups(
tumor,
groupby="condition",
groups=["BRCA1"],
reference="BRCA2",
method="wilcoxon"
)

deg = sc.get.rank_genes_groups_df(
tumor,
group="BRCA1"
)

deg.to_excel(
f"{BASE}/tables/BRCA1_vs_BRCA2_DEG.xlsx",
index=False
)

print("DEG table saved")

############################################################
# VOLCANO PLOT
############################################################

plt.figure(figsize=(7,6))

plt.scatter(
deg["logfoldchanges"],
-np.log10(deg["pvals_adj"] + 1e-300),
s=8
)

plt.xlabel("log2 Fold Change")
plt.ylabel("-log10 adjusted p-value")

plt.title("BRCA1 vs BRCA2 Differential Expression")

plt.axvline(1, linestyle="--")
plt.axvline(-1, linestyle="--")

plt.tight_layout()

plt.savefig(
f"{BASE}/plots/BRCA1_vs_BRCA2_volcano.png"
)

plt.close()

############################################################
# AUTOMATIC IMMUNE CELL DETECTION
############################################################

print("\nAnalyzing immune microenvironment")

immune_keywords = [
"cd","nk","macro","mono","dc","lymph","plasma"
]

immune_types = [
x for x in adata.obs["majority_voting"].unique()
if any(k in x.lower() for k in immune_keywords)
]

print("Detected immune cell types:")
print(immune_types)

immune = adata[
adata.obs["majority_voting"].isin(immune_types)
].copy()

immune_counts = pd.crosstab(
immune.obs["majority_voting"],
immune.obs["condition"]
)

immune_props = immune_counts.div(
immune_counts.sum(axis=0),
axis=1
) * 100

immune_props.to_excel(
f"{BASE}/tables/immune_cell_composition.xlsx"
)

############################################################
# IMMUNE MICROENVIRONMENT PLOT
############################################################

plt.figure(figsize=(12,7))

immune_props.T.plot(
kind="bar",
stacked=True
)

plt.ylabel("Immune cell proportion (%)")
plt.title("Tumor Immune Microenvironment")

plt.legend(
bbox_to_anchor=(1.02,1),
loc="upper left"
)

plt.tight_layout()

plt.savefig(
f"{BASE}/plots/immune_microenvironment.png"
)

plt.close()

############################################################
# UMAP VISUALIZATION
############################################################

print("\nGenerating UMAP plots")

sc.settings.figdir = f"{BASE}/plots"

sc.pl.umap(
adata,
color="majority_voting",
save="_celltypes.png",
show=False
)

sc.pl.umap(
adata,
color="condition",
save="_condition.png",
show=False
)

############################################################

print("\nAll interpretation analyses completed successfully!")