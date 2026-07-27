expression_data <- read.table(
  "raw_data.csv",
  sep = ",",
  header = TRUE,
  check.names = FALSE,
  row.names = 1
)
sample_metadata <- read.table(
  "lasso_group.csv",
  sep = ",",
  header = TRUE,
  check.names = FALSE,
  row.names = 1
)
gene_filter_table <- read.table(
  "after_intersection_mito.csv",
  sep = ",",
  header = TRUE,
  check.names = FALSE
)

genes_to_keep <- gene_filter_table$gene
expression_data <- expression_data[
  rownames(expression_data) %in% genes_to_keep,
]
lasso_data <- as.data.frame(t(expression_data))
lasso_data$Group <- sample_metadata[rownames(lasso_data), "group"]
lasso_data <- lasso_data[,
  c("Group", setdiff(names(lasso_data), "Group"))
]

library(glmnet)

predictor_matrix <- as.matrix(lasso_data[1:15])
response_vector <- lasso_data[, 1]

lasso_model <- glmnet(
  predictor_matrix,
  response_vector,
  family = "binomial",
  alpha = 1
)
print(lasso_model)
plot(lasso_model, xvar = "lambda", label = FALSE, sign.lambda = 1)
plot(lasso_model, xvar = "lambda")

set.seed(123456)
cv_model <- cv.glmnet(
  predictor_matrix,
  response_vector,
  family = "binomial",
  alpha = 1,
  nfolds = 5
)
plot(cv_model, sign.lambda = 1)

lambda_min <- cv_model$lambda.min
lambda_min
lambda_1se <- cv_model$lambda.1se
lambda_1se

coefficient_matrix_1se <- coef(lasso_model, s = lambda_1se)
coefficient_matrix_1se
write.csv(
  as.matrix(coefficient_matrix_1se),
  "lassos_1se.csv"
)

coefficient_matrix_min <- coef(lasso_model, s = lambda_min)
coefficient_matrix_min
write.csv(
  as.matrix(coefficient_matrix_min),
  "lassos_min.csv"
)
write.csv(
  data.frame(lambda_1se = lambda_1se, lambda_min = lambda_min),
  "lassos_lamda.csv"
)

library(randomForest)
library(pROC)
library(caret)
library(limma)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(patchwork)

set.seed(123456)

expression_data <- read.table(
  "raw_data.csv",
  header = TRUE,
  sep = ",",
  check.names = FALSE,
  row.names = 1
)
sample_metadata <- read.table(
  "RF_design.csv",
  header = TRUE,
  sep = ",",
  check.names = FALSE,
  row.names = 1
)

if (!"type" %in% colnames(sample_metadata)) {
  stop("数据框 'sample_metadata' 中没有找到 'type' 列。")
}

if (!all(sample_metadata$type %in% c("Control", "Disease"))) {
  stop("'type' 列包含非预期的值，应为 'Control' 或 'Disease'")
}

group_order <- ifelse(sample_metadata$type == "Control", 1, 2)
sorted_sample_names <- rownames(sample_metadata)[order(group_order)]
expression_data <- expression_data[, sorted_sample_names, drop = FALSE]

target_genes <- read.table(
  "after_intersection_mito.csv",
  header = TRUE,
  sep = ",",
  check.names = FALSE,
  row.names = 1
)
common_genes <- intersect(
  rownames(expression_data)[apply(expression_data, 1, function(values) {
    all(values != 5.02e-5)
  })],
  rownames(target_genes)
)
expression_data <- expression_data[common_genes, , drop = FALSE]
training_matrix <- as.matrix(t(expression_data))
colnames(training_matrix) <- make.names(colnames(training_matrix))
outcome_labels <- as.factor(sample_metadata[sorted_sample_names, "type"])
train_control <- trainControl(
  method = "cv",
  number = 5,
  savePredictions = TRUE,
  classProbs = TRUE,
  search = "grid"
)
tuning_grid <- expand.grid(
  mtry = floor(sqrt(ncol(training_matrix)))
)

set.seed(123456)
cv_model <- train(
  x = training_matrix,
  y = outcome_labels,
  method = "rf",
  ntree = 500,
  tuneGrid = tuning_grid,
  trControl = train_control
)

optimal_tree_count <- cv_model$finalModel$ntree
cat("最佳树的数量:", optimal_tree_count, "\n")

model_data <- as.data.frame(training_matrix)
model_data$outcome_label <- outcome_labels

rf_model <- randomForest(
  outcome_label ~ .,
  data = model_data,
  ntree = optimal_tree_count,
  importance = TRUE
)

pdf("最佳树的随机森林误差曲线.pdf", width = 8, height = 6)
plot(
  rf_model,
  main = paste("Random Forest with", optimal_tree_count, "Trees"),
  lwd = 2
)
dev.off()

predicted_classes <- predict(rf_model)
confusion_matrix <- confusionMatrix(predicted_classes, outcome_labels)

confusion_data <- as.data.frame(confusion_matrix$table)
confusion_data$Actual <- factor(
  confusion_data$Reference,
  levels = rev(levels(confusion_data$Reference))
)

confusion_plot <- ggplot(
  confusion_data,
  aes(x = Prediction, y = Actual, fill = Freq)
) +
  geom_tile(color = "white", alpha = 0.8) +
  geom_text(
    aes(label = Freq),
    color = "white",
    size = 6,
    fontface = "bold"
  ) +
  scale_fill_gradient(low = "#2196F3", high = "#F44336") +
  labs(
    title = paste(
      "Confusion Matrix\nAccuracy =",
      round(confusion_matrix$overall["Accuracy"], 3)
    ),
    x = "Predicted",
    y = "Actual"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  "RF_Confusion_Matrix.pdf",
  confusion_plot,
  width = 5,
  height = 4.5
)

predicted_probabilities <- predict(rf_model, type = "prob")[, 2]
model_roc <- roc(
  response = outcome_labels,
  predictor = predicted_probabilities
)

pdf("RF_Model_ROC_Curve.pdf", width = 6, height = 6)
plot(
  model_roc,
  main = paste("ROC Curve (AUC =", round(auc(model_roc), 3), ")"),
  col = "#1c61b6",
  lwd = 3,
  legacy.axes = TRUE
)
dev.off()

gene_importance <- data.frame(
  Gene = rownames(importance(rf_model)),
  Importance = importance(rf_model)[, "MeanDecreaseGini"],
  stringsAsFactors = FALSE
) %>%
  arrange(desc(Importance))

write.csv(
  gene_importance,
  "randomforest_Gene_Importance.csv",
  row.names = FALSE
)

top_genes <- head(gene_importance, 10)
write.csv(
  top_genes,
  "randomforest_最重要的10个与疾病相关的基因.csv",
  row.names = FALSE
)

importance_plot <- ggplot(
  top_genes,
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
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    axis.text.y = element_text(size = 11, face = "italic"),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10),
    panel.border = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

importance_plot

ggsave(
  "RF_Top10_Gene_Importance.pdf",
  importance_plot,
  width = 4,
  height = 4,
  device = "pdf"
)

gene_roc_data <- as.data.frame(training_matrix)
gene_roc_data$outcome_label <- outcome_labels

gene_roc_plots <- list()
for (gene_name in top_genes$Gene) {
  gene_roc <- roc(
    gene_roc_data$outcome_label,
    gene_roc_data[[gene_name]]
  )

  gene_roc_plots[[gene_name]] <- ggroc(
    gene_roc,
    legacy.axes = TRUE,
    color = "#1f78b4"
  ) +
    ggtitle(
      paste0(gene_name, "\nAUC = ", round(auc(gene_roc), 3))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5),
      panel.grid = element_blank()
    )
}

combined_roc_plot <- wrap_plots(gene_roc_plots, ncol = 5) +
  plot_annotation(
    title = "ROC Curves of Top 10 Genes",
    theme = theme(plot.title = element_text(hjust = 0.5))
  )

ggsave(
  "RF_Top10_Gene_ROCs.pdf",
  combined_roc_plot,
  width = 12,
  height = 6
)

sink("Model_Summary.txt")
cat("=== Random Forest Model Summary ===\n\n")
cat("Optimal number of trees:", optimal_tree_count, "\n\n")
cat("=== Cross-Validation Results ===\n")
print(cv_model$results)
cat("\n\n=== Final Model Performance ===\n")
print(confusion_matrix)
cat("\nModel AUC:", round(auc(model_roc), 3), "\n")
cat("\n=== Top 10 Important Genes ===\n")
print(top_genes)
sink()

