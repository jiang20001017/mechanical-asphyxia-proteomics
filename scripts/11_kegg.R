library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

differential_genes <- read.csv("input/genes_deg.csv")

gene_symbols <- differential_genes[[1]]
gene_symbols <- unique(gene_symbols[!is.na(gene_symbols) & gene_symbols != ""])

entrez_genes <- bitr(
    gene_symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Hs.eg.db"
)
cat("转换成功的基因数：", nrow(entrez_genes), "\n")

kegg_result <- enrichKEGG(
    entrez_genes$ENTREZID,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 0.05
)

kegg_result <- setReadable(
    kegg_result,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID"
)

kegg_table <- as.data.frame(kegg_result)
cat("KEGG 通路数：", nrow(kegg_table), "\n")

significant_kegg <- kegg_table %>%
    dplyr::filter(qvalue < 0.05) %>%
    dplyr::arrange(qvalue) %>%
    dplyr::mutate(
        GeneRatio = sapply(strsplit(as.character(GeneRatio), "/"), function(s) {
            as.numeric(s[1]) / as.numeric(s[2])
        }),
        pvalue = format(pvalue, scientific = FALSE),
        p.adjust = format(p.adjust, scientific = FALSE),
        qvalue = format(qvalue, scientific = FALSE)
    ) %>%
    dplyr::select(
        ID,
        Description,
        GeneRatio,
        BgRatio,
        pvalue,
        p.adjust,
        qvalue,
        geneID,
        Count
    )

write.csv(
    significant_kegg,
    "output/kegg_significant_pathways.csv",
    row.names = FALSE
)
cat("显著通路数：", nrow(significant_kegg), "\n")
