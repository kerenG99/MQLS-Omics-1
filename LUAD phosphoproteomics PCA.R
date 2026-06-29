library(reticulate)
library(ggplot2)

cptac <- import("cptac")
luad  <- cptac$Luad()

mat  <- py_to_r(luad$get_phosphoproteomics(source = "bcm"))
clin <- py_to_r(luad$get_clinical(source = "mssm"))

common <- intersect(rownames(mat), rownames(clin))
mat    <- mat[common, ]
stage  <- clin[common, "tumor_stage_pathological"]

keep  <- !is.na(stage)
mat   <- mat[keep, ]
stage <- stage[keep]
table(stage)

na_frac  <- apply(mat, 2, function(x) mean(is.na(x)))
mat      <- mat[, na_frac < 0.5, drop = FALSE]
mat      <- apply(mat, 2, function(x) { x[is.na(x)] <- median(x, na.rm = TRUE); x })
mat      <- as.matrix(mat)  # ensure matrix type after apply

# Remove zero/constant columns AFTER imputation (apply can reset variance)
col_vars <- apply(mat, 2, var)
mat      <- mat[, col_vars > 0, drop = FALSE]

# Extra safety: remove any remaining NA columns
mat <- mat[, colSums(is.na(mat)) == 0, drop = FALSE]

pca    <- prcomp(mat, center = TRUE, scale. = TRUE)
pca_df <- data.frame(pca$x[, 1:2], Stage = stage)

ggplot(pca_df, aes(PC1, PC2, color = Stage)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "Phosphoproteomics BCM - PCA by Tumor Stage") +
  theme_classic()

ggplot(pca_df, aes(PC1, PC2, color = Stage)) +
  geom_point(size = 1.5, alpha = 0.6) +
  xlim(quantile(pca_df$PC1, 0.01), quantile(pca_df$PC1, 0.99)) +
  ylim(quantile(pca_df$PC2, 0.01), quantile(pca_df$PC2, 0.99)) +
  labs(title = "Phosphoproteomics BCM - PCA Trimmed") +
  theme_classic()

ggplot(pca_df, aes(PC1, PC2, color = Stage)) +
  geom_point(size = 1.5, alpha = 0.6) +
  xlim(quantile(pca_df$PC1, 0.01), quantile(pca_df$PC1, 0.99)) +
  ylim(quantile(pca_df$PC2, 0.01), quantile(pca_df$PC2, 0.99)) +
  facet_wrap(~Stage) +
  labs(title = "Phosphoproteomics BCM - PCA Faceted by Stage") +
  theme_classic()
