SSRI Pharmacogenomics Modeling
Translational Transcriptomic Response Prediction (GSE146446)

Project Overview

This project develops predictive models of antidepressant treatment response using longitudinal transcriptomic data from a clinical study (GSE146446).

We investigate whether gene expression patterns can predict treatment response in patients treated with:

DLX (Duloxetine) – active antidepressant

PLB (Placebo) – control

The primary goal is to:

Identify response-associated transcriptional signatures

Compare baseline vs treatment-induced signals

Build predictive models for clinical response


Scientific Objective

Can gene expression predict antidepressant treatment response?

We test two hypotheses:

Baseline gene expression (T0) contains predictive signal

Treatment-induced transcriptional change (T8 − T0) improves prediction


📊 Dataset

GEO accession: GSE146446

Platform: Microarray transcriptomics

Total samples: 406

DLX-treated subjects: 96 (paired T0 and T8)

Outcome: Clinical response status (Responder vs Nonresponder)



Project Structure

ssri-pgx-mdd-gse146446/
│
├── data/
│   ├── raw/
│   └── processed/
│
├── scripts/
│   ├── 01_cleaning.R
│   ├── 02_DE_analysis.R
│   └── 03_modeling.R
│
├── results/
│   ├── DE_baseline_all.csv
│   ├── DE_delta_all.csv
│   ├── volcano_delta.png
│   ├── heatmap_delta.png
│   ├── ROC_baseline.png
│   ├── ROC_delta.png
│   └── model_*_metrics.txt
│
└── README.md


⚙️ How to Run

Step 1 — Set Working Directory
setwd("path_to/ssri-pgx-mdd-gse146446")
Step 2 — Data Processing & Dataset Construction
source("scripts/01_cleaning.R")

This script:

Downloads or loads GSE146446

Aligns expression matrix with clinical metadata

Cleans and standardizes variables

Builds:

Baseline modeling dataset (DLX, T0)

Delta dataset (DLX, T8 − T0 paired)

Outputs:

data/DLX_baseline_modeling_ready.rds

data/processed/DLX_delta_data.rds

Step 3 — Differential Expression Analysis
source("scripts/02_DE_analysis.R")

Performs:

Baseline DE (Responder vs Nonresponder)

Delta DE (T8 − T0 difference)

limma modeling

Multiple testing correction (BH FDR)

Outputs:

DE tables (.csv)

Volcano plot

Heatmap of top delta genes

Step 4 — Predictive Modeling
source("scripts/03_modeling.R")

Builds two models:

Model 1 — Baseline Expression

Input: T0 DLX samples

Feature selection: top variance genes (train-only)

Logistic regression

70/30 train-test split

ROC/AUC evaluation

Model 2 — Delta Expression

Input: paired T8 − T0 differences

Same modeling framework

ROC/AUC evaluation

Outputs:

ROC curves

Model metrics files



Methods

1️⃣ Clinical Data Harmonization

Standardized clinical fields:

visit

response

treatment

subject_id

Ensured sample alignment

Removed duplicates

Paired longitudinal samples

2️⃣ Dimensionality Reduction

Variance filtering

Top 5000 genes retained

Feature selection performed on training set only (to avoid leakage)

3️⃣ Differential Expression (limma)

Linear modeling

Responder vs Nonresponder contrast

Benjamini-Hochberg FDR correction

Visualization: volcano + heatmap

4️⃣ Predictive Modeling

Logistic regression

Train/test split (70/30)

ROC analysis using pROC

AUC calculation

Confusion matrix reporting


📈 Results

Model	AUC
Baseline	~0.57
Delta	~0.62
Key Finding

Treatment-induced transcriptional changes (delta) outperform baseline expression in predicting treatment response.

This suggests that pharmacodynamic signal provides stronger predictive information than static baseline expression.



Biological Interpretation

Minimal baseline differences between responders and nonresponders

Stronger separation observed in treatment-induced changes

Supports the hypothesis that dynamic response reflects drug mechanism



Reproducibility

The entire pipeline can be reproduced with:

source("scripts/01_cleaning.R")
source("scripts/02_DE_analysis.R")
source("scripts/03_modeling.R")

No manual steps required.



Skills Demonstrated

GEO data retrieval (GEOquery)

Clinical metadata harmonization

Longitudinal delta modeling

Differential expression (limma)

Feature selection without leakage

Logistic regression modeling

ROC/AUC evaluation

Reproducible project structuring


Future Improvements

5-fold cross-validation

Random Forest comparison

Pathway enrichment analysis

External validation dataset


Industry Relevance

This project demonstrates competencies relevant to:

Bioinformatics Associate

Translational Data Analyst

Biomarker Discovery Scientist

Clinical analytics roles in pharma or CRO


Author

**** Eugenia Yi ****

Translational Pharmacogenomics Modeling Project
Built in R
