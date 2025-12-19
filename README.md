# Metabolic-Heterogeneity-in-MTC
This project presents metabolomic (LC/GC-MS) and transcriptomic (RNA-seq) data from a cohort of 101 medullary thyroid cancer (MTC) patients. The primary goals were to define metabolic subtypes in MTC and identify clinically relevant biomarkers.
[README1.md]

> File names can be adjusted as needed.  
> Raw data are **not included** in this repository.

---

## Main Analyses

### 1. GSVA-based Metabolic Pathway Scoring

- KEGG metabolic pathways are used as gene sets
- GSVA is applied to normalized gene expression data
- Pathways are filtered based on variability (standard deviation)
- Z-score normalization is performed before clustering

**Key output**
- GSVA score matrix (`gsva_data.txt`)
- Metabolic gene set object (`cellmarker_ssGSEA.Rdata`)

---

### 2. Consensus Clustering

- Unsupervised clustering based on GSVA scores
- Algorithm: k-means
- Distance metric: Pearson correlation
- Multiple cluster numbers (K = 2–8) evaluated
- Consensus matrices and clustering stability visualized

**Key output**
- Consensus clustering PDFs
- Sample group assignments (`group.Rda`)

---

### 3. Metabolomics Analysis

- Principal Component Analysis (PCA)
- Orthogonal Partial Least Squares Discriminant Analysis (OPLS-DA)
- 2D and 3D score visualization

This step is used to evaluate metabolic differences between predefined sample groups.

---

### 4. LASSO-based Feature Selection

- Logistic LASSO regression using `glmnet`
- 10-fold cross-validation
- Identification of features with non-zero coefficients
- Suitable for downstream modeling or biomarker screening

**Key output**
- Selected feature list (`lassoGene` object)

---

## Input Data Format

### Transcriptomics

- **Gene expression matrix**
  - Rows: genes
  - Columns: samples
  - Normalized values (e.g., FPKM)

- **Raw count matrix**
  - Used for supplementary analyses if needed

- **Metabolic gene sets**
  - CSV format
  - Columns:
    - Pathway / Metagene
    - Category / Cell type

# Multi-omics Metabolism Analysis Pipeline

This repository contains an R-based analysis pipeline for metabolism-focused multi-omics data analysis, including transcriptomics-derived pathway activity scoring, unsupervised clustering, metabolomics analysis, and feature selection using LASSO regression.

The code is designed for exploratory and downstream analysis in cancer metabolism studies and can be adapted to other disease contexts.

---

## Repository Structure


---

### Metabolomics

- **Metabolite abundance matrix**
  - Rows: metabolites
  - Columns: samples

- **Sample annotation file**
  - Sample IDs
  - Group labels

---

## Requirements

### R version
- R ≥ 4.0 recommended

### R packages

r
tidyverse
data.table
GSVA
ConsensusClusterPlus
sweep
ropls
glmnet
ggplot2
ggrepel
ggthemes



[README2.md]
# ResNet-based Deep Learning Model

This repository contains a Jupyter Notebook implementing a ResNet-based deep learning model for supervised learning tasks.

The notebook covers data loading, model construction, training, and evaluation in an end-to-end workflow.

---

## File Description


---

## Overview

The notebook includes the following steps:

1. Data loading and preprocessing  
2. ResNet model definition  
3. Model training  
4. Model evaluation and result visualization  

The code is intended for research and experimental purposes and can be adapted to different datasets.

---

## Requirements

### Environment

- Python ≥ 3.8
- Jupyter Notebook / JupyterLab

### Python Packages

Typical dependencies include:

```text
numpy
pandas
torch / tensorflow
matplotlib
scikit-learn
