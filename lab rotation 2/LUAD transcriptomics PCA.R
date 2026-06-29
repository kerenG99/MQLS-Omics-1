#TASK 2
#dowload cpatc
#dowload the cancer lung adenocarcinoma
#each has different omics data, transciptome , proteomics
#each omics data has stages- PCA plot for this
#cluster separatiion for each omics data
#filter

library(reticulate)
#accessing the cptac package
py_require("cptac")
library(ggplot2)

#creating cptac to hold the data
cptac <- import("cptac")
#creating luad dataset from cptac 
luad  <- cptac$Luad()
luad$list_data_sources()

# Get trans data from specific lab ie source is washington university,converting the python objet to an r object
transcriptomics <- py_to_r(luad$get_transcriptomics(source = "washu"))
#view
View(py_to_r(luad$get_transcriptomics(source = "washu")))

#extract clinical data
clin  <- py_to_r(luad$get_clinical(source = "mssm"))
#view
View(clin)

# storing data and creating objects to use
trans <- py_to_r(luad$get_transcriptomics(source = "washu"))
clin  <- py_to_r(luad$get_clinical(source = "mssm"))

# stage column
grep("stage", colnames(clin), value = TRUE, ignore.case = TRUE)

# Common patients in both datasets, keeping all columns in trans matrix and only stage column in clin
common <- intersect(rownames(trans), rownames(clin))
mat    <- trans[common, ]
stage  <- clin[common, "tumor_stage_pathological"]
View(mat)

# removing NA stages
keep  <- !is.na(stage)
mat   <- mat[keep, ]
stage <- stage[keep]
table(stage) 

# removing genes with no variation in mat dataset columns since rows is patients, and ignoring missing values
#greating a gene vars object as a result in order to keep genes with a variance greater than zero
gene_vars <- apply(mat, 2, var, na.rm = TRUE)
mat <- mat[, gene_vars > 0]
View(mat)

# PCA
pca    <- prcomp(mat, center = TRUE, scale. = TRUE)
pca_df <- data.frame(pca$x[, 1:2], Stage = stage)
View(pca_df) 

ggplot(pca_df, aes(PC1, PC2, color = Stage)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_classic()

# remove outliers
ggplot(pca_df, aes(PC1, PC2, color = Stage)) +
  geom_point(size = 1.5, alpha = 0.6) +
  xlim(-150, 100) +
  ylim(-100, 150) +
  labs(title = "LUAD Transcriptomics PCA by Tumor Stage") +
  theme_classic()

ggplot(pca_df, aes(PC1, PC2, color = Stage)) +
  geom_point(size = 1.5, alpha = 0.6) +
  xlim(-150, 100) +
  ylim(-100, 150) +
  facet_wrap(~Stage) + #miniplot for each stage
  labs(title = "LUAD Transcriptomics PCA by Tumor Stage") +
  theme_classic()
