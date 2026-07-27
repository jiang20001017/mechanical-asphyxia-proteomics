# Analysis Code for Mechanical Asphyxia Proteomics Study

This repository contains the analysis code used in the study:

**Combined Machine Learning Identification and Experimental Validation of NDUFS8 and SUCLG1 as Potential Biomarkers for Mechanical Asphyxia**

## Contents

- `LASSO_RF.R` — machine-learning analysis: LASSO logistic regression and random forest, used for candidate-protein prioritization.

Additional scripts for proteomic data processing, differential-expression analysis, bioinformatic analysis, and figure generation will be added.

## Data availability

The mass-spectrometry data have been deposited in the ProteomeXchange Consortium via the iProX partner repository under accession number **PXD081244** (iProX identifier: **IPX0018430000**).

## Reproducibility

- R version: 4.4.1
- Random seed: set.seed(123456)
- Key R packages:
  - glmnet 4.1-10
  - randomForest 4.7-1.2
  - caret 7.0-1
