# Analysis Code for Mechanical Asphyxia Proteomics Study

This repository contains the analysis code used in the study:

**Combined Machine Learning Identification and Experimental Validation of NDUFS8 and SUCLG1 as Potential Biomarkers for Mechanical Asphyxia**

## Repository structure

- `scripts/` — all analysis scripts, numbered in the order they are run
- `input/` — input files required by the scripts

## Analysis scripts (run in order)

1. `01_filter_genes.R` — missing-value filtering of the protein-abundance matrix
2. `02_dma_bi_statistics.R` — differential-abundance analysis, MA vs BI (log2 scale)
3. `03_dma_hs_statistics.R` — differential-abundance analysis, MA vs HS (log2 scale)
4. `04_pca.R` — principal component analysis
5. `05_heatmap.R` — heatmap of the differentially abundant proteins
6. `06a_volcano_dma_bi.R` / `06b_volcano_dma_hs.R` — volcano plots
7. `07_lasso.R` — LASSO logistic regression
8. `08_elastic_net.R` — elastic-net regression
9. `09_random_forest.R` — random forest analysis
10. `10_go_ridgeplot.R` — GO enrichment analysis
11. `11_kegg.R` — KEGG pathway enrichment analysis
12. `12_gsea.R` — single-protein gene set enrichment analysis (GSEA)

## Data availability

The mass-spectrometry proteomics data have been deposited in the ProteomeXchange Consortium via the iProX partner repository under accession number **PXD081244** (iProX identifier: IPX0018430000). Owing to the confidentiality agreements and data-protection policy of our laboratory, the complete protein-abundance matrices will be released publicly upon publication; the processed input matrix required to run the machine-learning analysis is provided in `input/`.

## Reproducibility

- R version: 4.4.1
- Random seed: `set.seed(123456)`
- Key R packages:
  - glmnet 4.1-10
  - randomForest 4.7-1.2
  - caret 7.0-1
  - clusterProfiler, enrichplot, org.Hs.eg.db (GO/KEGG/GSEA)
  - ggplot2, pheatmap, ggbiplot (visualization)
