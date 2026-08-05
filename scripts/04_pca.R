library(ggplot2)
library(ggbiplot)

expression_matrix <- read.table(
    "output/expression_filtered.csv",
    header = TRUE,
    row.names = 1,
    sep = ",",
    check.names = FALSE
)

expression_matrix <- log2(expression_matrix)
expression_matrix <- t(expression_matrix)
expression_matrix <- expression_matrix[,
    apply(expression_matrix, 2, function(x) var(x) > 0),
    drop = FALSE
]

sample_groups <- read.csv("input/group_three.csv", header = TRUE, row.names = 1)
sample_groups <- sample_groups[rownames(expression_matrix), , drop = FALSE]

if (any(is.na(sample_groups[, 1]))) {
    stop("分组文件与样本名对不上，请检查 group_three.csv 的行名。")
}

pca_result <- prcomp(expression_matrix, scale. = TRUE)

pca_plot <- ggbiplot(
    pca_result,
    var.axes = FALSE,
    obs.scale = 1,
    groups = sample_groups[, 1],
    ellipse = TRUE,
    circle = FALSE
) +
    geom_text(
        aes(label = rownames(expression_matrix)),
        vjust = 1.5,
        size = 3.5
    ) +
    theme_bw() +
    ggtitle("Principal Component Analysis (PCA)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(
        values = c(DMA = "#AA0000", HS = "#005696", BI = "#0E9E57")
    ) +
    theme(
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank()
    )

ggsave("output/pca_plot.pdf", pca_plot, width = 7, height = 6)
