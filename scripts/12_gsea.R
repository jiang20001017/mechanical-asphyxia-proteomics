library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(ggplot2)
library(stringr)

expression_matrix <- read.csv(
    "output/expression_filtered.csv",
    row.names = 1,
    check.names = FALSE
)
expression_matrix <- as.matrix(expression_matrix)
expression_matrix <- log2(expression_matrix)

hallmark_sets <- read.gmt("input/hallmark.gmt")

run_single_gene_gsea <- function(target, expression_matrix, gmt, outprefix) {
    target_expression <- expression_matrix[target, ]

    gene_correlations <- apply(expression_matrix, 1, function(gene_expression) {
        suppressWarnings(cor(
            target_expression,
            gene_expression,
            method = "spearman"
        ))
    })
    gene_correlations <- gene_correlations[names(gene_correlations) != target]
    gene_correlations <- gene_correlations[is.finite(gene_correlations)]
    geneList <- sort(gene_correlations, decreasing = TRUE)

    set.seed(123)
    gsea_result <- GSEA(
        geneList,
        TERM2GENE = gmt,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        minGSSize = 10,
        maxGSSize = 2000,
        eps = 0
    )

    gsea_result@result$Description <- gsea_result@result$Description %>%
        str_remove("^HALLMARK_") %>%
        str_replace_all("_", " ")

    write.csv(
        gsea_result@result,
        paste0(outprefix, "_results.csv"),
        row.names = FALSE
    )

    if (nrow(gsea_result@result) == 0) {
        message(target, "：未富集到任何通路")
        return(gsea_result)
    }

    top5 <- gsea_result@result %>%
        arrange(p.adjust) %>%
        dplyr::slice(1:5)
    top5_id <- top5$ID

    new_lab <- paste0(
        top5$Description,
        "  (q = ",
        formatC(top5$p.adjust, format = "e", digits = 2),
        ")"
    )
    idx <- match(top5_id, gsea_result@result$ID)
    gsea_result@result$Description[idx] <- new_lab

    p_main <- gseaplot2(gsea_result, geneSetID = top5_id, pvalue_table = FALSE)
    ggsave(paste0(outprefix, "_top5.pdf"), p_main, width = 8, height = 6)

    return(gsea_result)
}

gsea_NDUFS8 <- run_single_gene_gsea(
    "NDUFS8",
    expression_matrix,
    hallmark_sets,
    "output/gsea_ndufs8"
)
gsea_SUCLG1 <- run_single_gene_gsea(
    "SUCLG1",
    expression_matrix,
    hallmark_sets,
    "output/gsea_suclg1"
)

head(gsea_NDUFS8@result[, c("Description", "NES", "pvalue", "p.adjust")], 10)
head(gsea_SUCLG1@result[, c("Description", "NES", "pvalue", "p.adjust")], 10)
