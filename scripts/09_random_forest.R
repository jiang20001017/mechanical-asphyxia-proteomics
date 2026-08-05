library(randomForest)
library(pROC)
library(caret)
library(limma)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(patchwork)
set.seed(123456)
expression_matrix <- read.table(
    "input/expression_candidate.csv",
    header = T,
    sep = ",",
    check.names = F,
    row.names = 1
)
sample_design <- read.table(
    "input/group_rf.csv",
    header = T,
    sep = ",",
    check.names = F,
    row.names = 1
)
expression_matrix <- log2(expression_matrix)
if (!"type" %in% colnames(sample_design)) {
    stop("缺少 type 列")
}
if (!all(sample_design$type %in% c("Control", "Disease"))) {
    stop("type 值异常")
}
group_order <- ifelse(sample_design$type == "Control", 1, 2)
sorted_colnames <- rownames(sample_design)[order(group_order)]
expression_matrix <- expression_matrix[, sorted_colnames, drop = FALSE]
expression_matrix <- expression_matrix[
    sort(rownames(expression_matrix)),
    ,
    drop = FALSE
]
expression_matrix <- as.matrix(t(expression_matrix))
colnames(expression_matrix) <- make.names(colnames(expression_matrix))
outcome <- as.factor(sample_design[sorted_colnames, "type"])
cv_control <- trainControl(
    method = "cv",
    number = 5,
    savePredictions = TRUE,
    classProbs = TRUE,
    search = "grid"
)
tuning_grid <- expand.grid(mtry = c(floor(sqrt(ncol(expression_matrix)))))
set.seed(123456)
cv_model <- train(
    x = expression_matrix,
    y = outcome,
    method = "rf",
    ntree = 500,
    tuneGrid = tuning_grid,
    trControl = cv_control
)
best_ntree <- cv_model$finalModel$ntree
cat("最佳树的数量:", best_ntree, "\n")
final_model <- randomForest(
    outcome ~ .,
    data = expression_matrix,
    ntree = 500,
    importance = TRUE,
    nodesize = 3
)
pdf("output/random_forest_error_curve.pdf", width = 8, height = 6)
plot(final_model, main = "Random Forest with 500 Trees", lwd = 2)
dev.off()
importance_df <- data.frame(
    Gene = rownames(importance(final_model)),
    Importance = importance(final_model)[, "MeanDecreaseGini"],
    stringsAsFactors = FALSE
) %>%
    arrange(desc(Importance))
write.csv(
    importance_df,
    "output/random_forest_gene_importance.csv",
    row.names = FALSE
)
top10 <- head(importance_df, 10)
write.csv(top10, "output/random_forest_top10_genes.csv", row.names = FALSE)
importance_plot <- ggplot(
    top10,
    aes(x = reorder(Gene, Importance), y = Importance, fill = Importance)
) +
    geom_bar(stat = "identity", width = 0.7) +
    coord_flip() +
    scale_fill_gradient(low = "#AFCBE3", high = "#F2B1B2") +
    labs(
        x = "Gene",
        y = "Importance (MeanDecreaseGini)",
        title = "Top 10 Important Genes"
    ) +
    theme_bw() +
    theme(
        axis.text.y = element_text(size = 11, face = "italic"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
ggsave(
    "output/random_forest_top10_importance.pdf",
    importance_plot,
    width = 4,
    height = 4
)
