############################################################
# COMPLETE SINGLE-CELL RNASEQ PIPELINE
############################################################

import os
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scrublet as scr
import gseapy as gp
import celltypist
import warnings

warnings.filterwarnings("ignore")

############################################################
# 1. CREATE OUTPUT DIRECTORIES
############################################################

BASE_DIR = "scRNA_results"

dirs = [
    f"{BASE_DIR}/data",
    f"{BASE_DIR}/plots/QC",
    f"{BASE_DIR}/plots/UMAP",
    f"{BASE_DIR}/plots/DEG",
    f"{BASE_DIR}/plots/pathways",
    f"{BASE_DIR}/plots/celltypes",
    f"{BASE_DIR}/tables/DEG",
    f"{BASE_DIR}/tables/pathways",
]

for d in dirs:
    os.makedirs(d, exist_ok=True)

print("Directories created")

############################################################
# 2. LOAD ALL H5 FILES
############################################################

data_dir = "."
adatas = []

for file in os.listdir(data_dir):

    if file.endswith(".h5"):

        print("Loading:", file)

        adata = sc.read_10x_h5(file)

        adata.var_names_make_unique()

        sample = file.replace(".h5","")
        name = sample.upper()

        if "BRCA1" in name:
            condition = "BRCA1"
        elif "BRCA2" in name:
            condition = "BRCA2"
        elif "NORMAL" in name:
            condition = "Normal"
        else:
            condition = "Unknown"

        if "GSM502" in name:
            batch = "Normal_batch"
        elif "GSM699" in name:
            batch = "Tumor_batch"
        else:
            batch = "Unknown"

        adata.obs["sample"] = sample
        adata.obs["condition"] = condition
        adata.obs["batch"] = batch

        adatas.append(adata)

############################################################
# 3. MERGE DATASETS
############################################################

adata = ad.concat(adatas, join="outer")

adata.obs_names_make_unique()
adata.var_names_make_unique()

print("Combined dataset:", adata)

adata.write_h5ad(f"{BASE_DIR}/data/BC_data_combined.h5ad")

############################################################
# 4. QC METRICS
############################################################

adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=[20],
    inplace=True
)

sc.pl.violin(
    adata,
    ["n_genes_by_counts","total_counts","pct_counts_mt"],
    show=False
)

plt.savefig(f"{BASE_DIR}/plots/QC/QC_violin.png")
plt.close()

############################################################
# 5. BASIC FILTERING
############################################################

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

############################################################
# 6. DOUBLETS
############################################################

print("Running doublet detection")

X_counts = adata.X

if not isinstance(X_counts,np.ndarray):
    X_counts = X_counts.toarray()

scrub = scr.Scrublet(X_counts)

scores, preds = scrub.scrub_doublets()

adata.obs["doublet_score"] = scores
adata.obs["doublet"] = preds

adata = adata[~adata.obs["doublet"]].copy()

print("Cells after filtering:", adata.n_obs)

adata.write_h5ad(f"{BASE_DIR}/data/BC_data_preprocessed.h5ad")

############################################################
# 7. NORMALIZATION
############################################################

adata.raw = adata

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    subset=True
)

sc.pp.scale(adata,max_value=10)

############################################################
# 8. PCA + UMAP
############################################################

sc.tl.pca(adata)

sc.pp.neighbors(adata,n_neighbors=15,n_pcs=30)

sc.tl.umap(adata)

sc.tl.leiden(adata,resolution=0.5)

sc.pl.umap(
    adata,
    color=["leiden","condition"],
    show=False
)

plt.savefig(f"{BASE_DIR}/plots/UMAP/UMAP_clusters.png")
plt.close()

adata.write_h5ad(f"{BASE_DIR}/data/BC_data_analyzed.h5ad")

############################################################
# 9. DIFFERENTIAL EXPRESSION
############################################################

print("Running DEG analysis")

sc.tl.rank_genes_groups(
    adata,
    "condition",
    method="wilcoxon"
)

deg = sc.get.rank_genes_groups_df(adata,group=None)

deg.to_excel(
    f"{BASE_DIR}/tables/DEG/condition_DEG.xlsx",
    index=False
)

############################################################
# 10. VOLCANO PLOT
############################################################

plt.figure(figsize=(6,5))

plt.scatter(
    deg["logfoldchanges"],
    -np.log10(deg["pvals_adj"]),
    s=6
)

plt.xlabel("log2 Fold Change")
plt.ylabel("-log10 adj p-value")

plt.savefig(f"{BASE_DIR}/plots/DEG/volcano.png")
plt.close()

############################################################
# 11. PATHWAY ENRICHMENT
############################################################

print("Running pathway enrichment")

ranking = deg[["names","logfoldchanges"]]
ranking.columns=["gene","score"]

ranking = ranking.sort_values("score",ascending=False)

pre_res = gp.prerank(
    rnk=ranking.set_index("gene")["score"],
    gene_sets="KEGG_2021_Human",
    outdir=None
)

pathways = pre_res.res2d

pathways.to_excel(
    f"{BASE_DIR}/tables/pathways/kegg_enrichment.xlsx",
    index=False
)

############################################################
# 12. PATHWAY PLOT
############################################################

top = pathways.head(10)

plt.figure(figsize=(6,4))

plt.barh(top["Term"],top["NES"])

plt.xlabel("NES")

plt.tight_layout()

plt.savefig(f"{BASE_DIR}/plots/pathways/top_pathways.png")
plt.close()

############################################################
# 13. CELLTYPE ANNOTATION
############################################################

print("Running CellTypist annotation")

adata.obs_names_make_unique()

pred = celltypist.annotate(
    adata,
    model="Cells_Adult_Breast.pkl",
    majority_voting=True
)

adata = pred.to_adata()

sc.pl.umap(
    adata,
    color="majority_voting",
    show=False
)

plt.savefig(f"{BASE_DIR}/plots/celltypes/celltypes_umap.png")
plt.close()

adata.write_h5ad(
    f"{BASE_DIR}/data/BC_data_celltypist.h5ad"
)

############################################################

print("Pipeline completed successfully")