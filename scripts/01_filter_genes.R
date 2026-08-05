input_file <- "input/expression_all.csv"
output_file <- "output/expression_filtered.csv"
fill_value <- 5.02e-05
max_fill_per_group <- 1L
dat <- read.csv(input_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
gene_col <- "sample"
dma_cols <- grep("^DMA[0-9]+$", names(dat), value = TRUE)
bi_cols <- grep("^BI[0-9]+$", names(dat), value = TRUE)
hs_cols <- grep("^HS[0-9]+$", names(dat), value = TRUE)
if (!gene_col %in% names(dat)) {
    stop("找不到基因列：", gene_col)
}
if (length(dma_cols) == 0L || length(bi_cols) == 0L || length(hs_cols) == 0L) {
    stop("未找到 DMA、BI 或 HS 样本列，请检查列名。")
}
value_cols <- c(dma_cols, bi_cols, hs_cols)
original_dat <- dat
dat[value_cols] <- lapply(dat[value_cols], function(x) {
    x_num <- suppressWarnings(as.numeric(x))
    if (any(is.na(x_num) & !is.na(x) & nzchar(trimws(x)))) {
        stop("表达量列中发现不能转换为数值的内容。")
    }
    x_num
})
is_fill <- function(x) {
    !is.na(x) & abs(x - fill_value) <= .Machine$double.eps^0.5 * max(1, abs(fill_value))
}
dma_fill_n <- rowSums(sapply(dat[dma_cols], is_fill))
bi_fill_n <- rowSums(sapply(dat[bi_cols], is_fill))
hs_fill_n <- rowSums(sapply(dat[hs_cols], is_fill))
keep <- dma_fill_n <= max_fill_per_group & bi_fill_n <= max_fill_per_group & hs_fill_n <= max_fill_per_group
result <- original_dat[keep, , drop = FALSE]
rownames(result) <- NULL
write.csv(result, output_file, row.names = FALSE, na = "")
message("筛选完成：", nrow(result), " 个基因")
message("结果文件：", output_file)
