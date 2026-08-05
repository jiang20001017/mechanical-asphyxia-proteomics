library(glmnet)

expression_matrix <- read.table(
    "input/expression_candidate.csv",
    sep = ",",
    header = TRUE,
    check.names = FALSE,
    row.names = 1
)
sample_groups <- read.table(
    "input/group_binary.csv",
    sep = ",",
    header = TRUE,
    check.names = FALSE,
    row.names = 1
)

expression_matrix <- log2(expression_matrix)
expression_matrix <- as.data.frame(t(expression_matrix))
expression_matrix$Group <- sample_groups[rownames(expression_matrix), "group"]
expression_matrix <- expression_matrix[,
    c("Group", setdiff(names(expression_matrix), "Group"))
]

predictor_matrix <- as.matrix(expression_matrix[, -1])
outcome <- expression_matrix[, 1]

lasso_model <- glmnet(predictor_matrix, outcome, family = "binomial", alpha = 1)
print(lasso_model)

set.seed(123456)
cv_model <- cv.glmnet(
    predictor_matrix,
    outcome,
    family = "binomial",
    alpha = 1,
    nfolds = 5
)

pdf("output/lasso_coefficient_path.pdf", width = 7, height = 6)
plot(lasso_model, xvar = "lambda")
dev.off()

pdf("output/lasso_cv_curve.pdf", width = 7, height = 6)
plot(cv_model)
dev.off()

lambda_min <- cv_model$lambda.min
lambda_1se <- cv_model$lambda.1se
print(lambda_min)
print(lambda_1se)

coefficients_1se <- coef(lasso_model, s = lambda_1se)
coefficients_min <- coef(lasso_model, s = lambda_min)
print(coefficients_1se)
print(coefficients_min)

write.csv(as.matrix(coefficients_1se), "output/lasso_coefficients_1se.csv")
write.csv(as.matrix(coefficients_min), "output/lasso_coefficients_min.csv")
write.csv(
    data.frame(lambda_1se = lambda_1se, lambda_min = lambda_min),
    "output/lasso_lambda_values.csv"
)
