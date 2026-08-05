input_file <- "output/expression_filtered.csv"
output_file <- "output/dma_hs_statistics.csv"
dat <- read.csv(
    input_file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8-BOM"
)
gene_col <- "sample"
dma_cols <- grep("^DMA[0-9]+$", names(dat), value = TRUE)
hs_cols <- grep("^HS[0-9]+$", names(dat), value = TRUE)
if (!gene_col %in% names(dat)) {
    stop("找不到基因列：", gene_col)
}
if (length(dma_cols) == 0L || length(hs_cols) == 0L) {
    stop("未找到 DMA 或 HS 样本列，请检查列名。")
}
value_cols <- c(dma_cols, hs_cols)
dat[value_cols] <- lapply(dat[value_cols], function(x) {
    x_num <- suppressWarnings(as.numeric(x))
    if (any(is.na(x_num) & !is.na(x) & nzchar(trimws(x)))) {
        stop("表达量列中发现不能转换为数值的内容。")
    }
    x_num
})
result <- dat
rownames(result) <- NULL
result_dma <- as.matrix(result[dma_cols])
result_hs <- as.matrix(result[hs_cols])
log2_dma <- log2(result_dma)
log2_hs <- log2(result_hs)
result$DMA_mean <- rowMeans(result_dma)
result$HS_mean <- rowMeans(result_hs)
result$logFC <- rowMeans(log2_dma) - rowMeans(log2_hs)
result$FC <- 2^result$logFC
calc_p <- function(i) {
    tryCatch(
        t.test(log2_dma[i, ], log2_hs[i, ], var.equal = TRUE)$p.value,
        error = function(e) NA_real_
    )
}
calc_d <- function(i) {
    x <- log2_dma[i, ]
    y <- log2_hs[i, ]
    s <- sqrt(
        ((length(x) - 1) * var(x) + (length(y) - 1) * var(y)) /
            (length(x) + length(y) - 2)
    )
    if (is.na(s) || s == 0) {
        return(NA_real_)
    }
    (mean(x) - mean(y)) / s
}
calc_ci <- function(i) {
    tryCatch(
        unname(t.test(log2_dma[i, ], log2_hs[i, ], var.equal = TRUE)$conf.int),
        error = function(e) c(NA_real_, NA_real_)
    )
}
result$pvalue <- vapply(seq_len(nrow(result)), calc_p, numeric(1))
result$Sig <- ifelse(
    !is.na(result$pvalue) & result$pvalue < 0.05 & result$logFC > 0.263,
    "Up",
    ifelse(
        !is.na(result$pvalue) &
            result$pvalue < 0.05 &
            result$logFC < -0.263,
        "Down",
        ""
    )
)
result$adj_pvalue_BH <- p.adjust(result$pvalue, method = "BH")
result$Cohen_d <- vapply(seq_len(nrow(result)), calc_d, numeric(1))
ci <- t(vapply(seq_len(nrow(result)), calc_ci, numeric(2)))
result$log2FC_CI95_low <- ci[, 1]
result$log2FC_CI95_high <- ci[, 2]
names(result)[names(result) == gene_col] <- "Gene_symbol"
result <- result[c(
    "Gene_symbol",
    dma_cols,
    hs_cols,
    "DMA_mean",
    "HS_mean",
    "FC",
    "pvalue",
    "logFC",
    "Sig",
    "adj_pvalue_BH",
    "Cohen_d",
    "log2FC_CI95_low",
    "log2FC_CI95_high"
)]
write.csv(result, output_file, row.names = FALSE, na = "")
message("统计完成：", nrow(result), " 个基因")
message("结果文件：", output_file)
