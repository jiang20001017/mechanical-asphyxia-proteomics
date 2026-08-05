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

set.seed(123456)
foldid <- sample(rep(1:5, length.out = nrow(predictor_matrix)))
alphas <- seq(0, 1, by = 0.1)

cv_models <- lapply(alphas, function(alpha_value) {
    cv.glmnet(
        predictor_matrix,
        outcome,
        family = "binomial",
        alpha = alpha_value,
        foldid = foldid
    )
})

best_alpha_index <- which.min(sapply(cv_models, function(cv) min(cv$cvm)))
best_alpha <- alphas[best_alpha_index]
cv_model <- cv_models[[best_alpha_index]]
print(best_alpha)

enet_model <- glmnet(
    predictor_matrix,
    outcome,
    family = "binomial",
    alpha = best_alpha
)

pdf("output/elastic_net_coefficient_path.pdf", width = 7, height = 6)
plot(enet_model, xvar = "lambda")
dev.off()

pdf("output/elastic_net_cv_curve.pdf", width = 7, height = 6)
plot(cv_model)
dev.off()

lambda_min <- cv_model$lambda.min
lambda_1se <- cv_model$lambda.1se
print(lambda_min)
print(lambda_1se)

coef_1se <- coef(enet_model, s = lambda_1se)
coef_min <- coef(enet_model, s = lambda_min)
print(coef_1se)
print(coef_min)

write.csv(as.matrix(coef_1se), "output/elastic_net_coefficients_1se.csv")
write.csv(as.matrix(coef_min), "output/elastic_net_coefficients_min.csv")
write.csv(
    data.frame(
        best_alpha = best_alpha,
        lambda_1se = lambda_1se,
        lambda_min = lambda_min
    ),
    "output/elastic_net_tuning_parameters.csv"
)
