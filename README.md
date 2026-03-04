# Single-Cell RNA-Seq Analysis of Breast Cancer Samples

Installation

Clone the repository:

git clone https://github.com/priteshparab/scRNA_Analysis.git
cd scRNA_Analysis

*Create a Python environment (recommended):*

conda create -n scrna_env python=3.10
conda activate scrna_env

# Install required dependencies:

pip install -r requirements.txt

Dataset source:

GEO Series: **GSE164898**

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164898

The dataset contains **single-cell transcriptomic profiles of human breast tissue samples**, including healthy controls and BRCA1/BRCA2 mutation carriers.

---

## Samples Used in This Study

### Normal Breast Tissue

| Sample ID | Description | Disease Status |
|-----------|-------------|---------------|
| GSM5022599 | Breast tissue | Healthy |
| GSM5022603 | Breast tissue | Healthy |

Experimental details:
- Organism: Homo sapiens
- Tissue: Breast
- Library strategy: RNA-Seq
- Library source: Transcriptomic
- Library selection: cDNA
- Sequencing platform: Illumina NovaSeq 6000

---

### BRCA2 Mutation Samples

| Sample ID | Tissue | Cell Type | Genotype |
|-----------|--------|-----------|----------|
| GSM6998337 | Breast | Epithelial | BRCA2 |
| GSM6998340 | Breast | Epithelial | BRCA2 |

Experimental details:
- Organism: Homo sapiens
- Tissue: Breast
- Library strategy: RNA-Seq
- Library source: Transcriptomic single cell
- Library selection: cDNA
- Sequencing platform: Illumina NovaSeq 6000

---

### BRCA1 Mutation Samples

| Sample ID | Tissue | Cell Type | Genotype |
|------------|--------|------------|-------|
| GSM6998341 | Breast | Epithelial | BRCA1 |
| GSM6998343 | Breast | Epithelial | BRCA1 |

Experimental details:
- Organism: Homo sapiens
- Tissue: Breast
- Library strategy: RNA-Seq
- Library source: Transcriptomic single cell
- Library selection: cDNA
- Sequencing platform: Illumina NovaSeq 6000

---


## Project Overview

This project performs **single-cell RNA sequencing (scRNA-seq) analysis** on breast cancer datasets containing:

* **Normal samples**
* **BRCA1 tumor samples**
* **BRCA2 tumor samples**

The goal of this analysis is to:

1. Identify **different cell types** present in the tissue.
2. Compare **tumor vs normal cell composition**.
3. Identify **genes that differ between conditions**.
4. Study **tumor microenvironment and immune infiltration**.
5. Perform **pathway enrichment analysis** to understand biological processes involved.

The analysis pipeline uses tools such as:

* **Scanpy** – core single-cell analysis
* **Scrublet** – doublet detection
* **CellTypist** – cell type annotation
* **GSEApy** – pathway enrichment analysis

---

# Folder Structure

The results are organized inside the `scRNA_results` directory.

```
scRNA_results
│
├── data
├── plots
├── tables
└── interpretation
```

Each folder contains specific outputs from different stages of the analysis.

---

# 1. DATA FOLDER

```
scRNA_results/data
```

This folder contains **processed single-cell datasets** saved in `.h5ad` format.

These files store the **gene expression matrix and metadata** used in each stage of the pipeline.

---

## BC_data_combined.h5ad

This file contains the **merged dataset of all samples**.

Information stored inside:

* Gene expression matrix
* Sample name
* Condition (Normal / BRCA1 / BRCA2)
* Batch information

At this stage:

* No filtering or normalization has been performed yet.

Purpose:

* Serves as the **starting point for downstream analysis**.

---

## BC_data_preprocessed.h5ad

This dataset contains **quality-controlled cells**.

Processing performed:

* Removal of **low-quality cells**
* Removal of **high mitochondrial content cells**
* **Doublet detection** using Scrublet
* Filtering genes expressed in very few cells

Purpose:

* Ensures only **high-quality single cells** remain for analysis.

---

## BC_data_analyzed.h5ad

This dataset contains the **core single-cell analysis results**.

Processing steps performed:

* Normalization
* Log transformation
* Highly variable gene selection
* PCA dimensionality reduction
* UMAP visualization
* Leiden clustering
* Marker gene identification

Purpose:

* Used to study **cell clusters and gene expression patterns**.

---

## BC_data_celltypist.h5ad

This dataset contains **cell type annotations**.

Cell types were assigned using **CellTypist**, which predicts cell identity based on known reference datasets.

Examples of predicted cell types:

* T cells
* B cells
* Macrophages
* Fibroblasts
* Epithelial cells

Purpose:

* Enables **biological interpretation of cell clusters**.

---

# 2. PLOTS FOLDER

```
scRNA_results/plots
```

This folder contains **visualizations generated during the analysis**.

---

## QC

Quality control plots showing:

* Number of genes per cell
* Total RNA counts
* Percentage of mitochondrial genes

Purpose:

* Helps identify and remove **low-quality cells**.

---

## UMAP

UMAP plots visualize **cell populations in 2D space**.

Cells that cluster together typically share similar gene expression patterns.

UMAP plots may show:

* Condition (Normal vs Tumor)
* Cluster identity
* Cell types

Purpose:

* Allows visual exploration of **cell heterogeneity**.

---

## DEG

Plots related to **Differential Gene Expression analysis**.

Includes:

* Volcano plots
* Gene ranking plots

Purpose:

* Shows which genes are **upregulated or downregulated** between conditions.

---

## pathways

Plots from **pathway enrichment analysis**.

These plots show biological pathways significantly enriched among differentially expressed genes.

Examples:

* Immune signaling pathways
* Cell cycle pathways
* DNA repair pathways

Purpose:

* Helps interpret **biological mechanisms driving tumor behavior**.

---

## celltypes

Plots related to **cell type composition and annotation**.

Includes:

* UMAP colored by predicted cell types
* Cell type distribution plots

Purpose:

* Understand **cellular composition of the tumor microenvironment**.

---

# 3. TABLES FOLDER

```
scRNA_results/tables
```

This folder contains **Excel files with numerical results**.

---

## DEG/condition_DEG.xlsx

Contains genes that are **differentially expressed between conditions**.

Typical comparisons include:

* BRCA1 vs Normal
* BRCA2 vs Normal
* BRCA1 vs BRCA2

Columns usually include:

* Gene name
* Log fold change
* Adjusted p-value

Purpose:

* Identify **genes associated with tumor progression**.

---

## pathways/kegg_enrichment.xlsx

Contains results of **KEGG pathway enrichment analysis**.

Each row corresponds to a biological pathway.

Information included:

* Pathway name
* Enrichment score
* Statistical significance

Purpose:

* Identify **biological processes affected in tumors**.

---

# 4. INTERPRETATION FOLDER

```
scRNA_results/interpretation
```

This folder contains **higher-level biological interpretation analyses**.

These analyses help understand **tumor biology and microenvironment changes**.

---

## interpretation/tables

### celltype_proportions.xlsx

Shows the **percentage of each cell type in each condition**.

Example:

| Cell Type   | Normal | BRCA1 | BRCA2 |
| ----------- | ------ | ----- | ----- |
| T cells     | 12%    | 25%   | 20%   |
| Fibroblasts | 10%    | 30%   | 28%   |

Purpose:

* Identify how **tumor composition changes compared to normal tissue**.

---

### immune_cell_composition.xlsx

Shows the distribution of **immune cell populations** across conditions.

Example cell types:

* T cells
* B cells
* Macrophages
* NK cells

Purpose:

* Study the **tumor immune microenvironment**.

---

### BRCA1_vs_BRCA2_DEG.xlsx

Contains genes that differ between **BRCA1 and BRCA2 tumor subtypes**.

Purpose:

* Identify **molecular differences between tumor subtypes**.

---

## interpretation/plots

### celltype_composition.png

Stacked bar plot showing **cell type distribution across conditions**.

Purpose:

* Visual comparison of **tumor vs normal cellular composition**.

---

### immune_microenvironment.png

Plot showing **immune cell infiltration across conditions**.

Purpose:

* Understand **immune response within tumors**.

---

### BRCA1_vs_BRCA2_volcano.png

Volcano plot showing genes differentially expressed between **BRCA1 and BRCA2 tumors**.

Genes far from the center are strongly differentially expressed.

---

### UMAP_celltypes.png

UMAP visualization with cells colored by predicted cell type.

Purpose:

* Visualize how **different cell populations cluster together**.

---

# Summary of Biological Insights

This analysis enables us to study:

1. **Cellular diversity** in breast cancer tissue.
2. Differences between **Normal, BRCA1, and BRCA2 samples**.
3. Changes in **tumor microenvironment and immune infiltration**.
4. Genes and pathways involved in **tumor progression**.

---

# Software Used

* Python
* Scanpy
* CellTypist
* Scrublet
* GSEApy
* Matplotlib
* Pandas

---

# Output Formats

The project generates several types of outputs:

| File Type | Purpose                                         |
| --------- | ----------------------------------------------- |
| `.h5ad`   | processed single-cell datasets                  |
| `.xlsx`   | tables of genes, pathways, and cell proportions |
| `.png`    | visualization figures                           |

---

# Notes

* All datasets are processed using **standard single-cell RNA-seq analysis workflows**.
* Results should be interpreted in the context of **biological validation and existing literature**.
* Additional analyses (e.g., cell-cell communication or CNV inference) can further extend this study.

---



