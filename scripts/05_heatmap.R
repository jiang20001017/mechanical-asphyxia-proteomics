library(pheatmap)
expr <- read.csv(
    "output/expression_filtered.csv",
    row.names = 1,
    check.names = FALSE
)
deg <- read.csv(
    "input/genes_deg.csv",
    check.names = FALSE
)
deg_genes <- deg[[1]]
deg_genes <- intersect(deg_genes, rownames(expr))
rt <- as.matrix(expr[deg_genes, , drop = FALSE])
cat("热图纳入基因数：", nrow(rt), "\n")
ann_col <- read.csv(
    "input/group_three.csv",
    row.names = 1
)
ann_col <- ann_col[colnames(rt), , drop = FALSE]
ann_colors <- list(group = c(DMA = "#F7C2CD", BI = "#FDD379", HS = "#A6DAEF"))
pheatmap(
    rt,
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    color = colorRampPalette(c("#005696", "white", "#AA0000"))(300),
    show_colnames = TRUE,
    show_rownames = FALSE,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    cutree_rows = 2,
    cutree_cols = 1,
    scale = "row",
    border_color = NA,
    fontsize_row = 4,
    treeheight_row = 10,
    treeheight_col = 6,
    cellwidth = 17,
    cellheight = 3,
    main = "Heatmap of Differentially Expressed Proteins",
    fontsize = 8,
    filename = "output/differential_expression_heatmap.pdf",
    width = 7,
    height = 12
)
