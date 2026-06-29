##installing packages

install.packages(c(
  "data.table",
  "ggplot2",
  "tidyverse",
  "reticulate"
))

#required biconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("MOFA2", "MOFAdata"))

#seting up python
library(reticulate)
install_miniconda()

#mofa backend
library(reticulate)

#create a fresh environment
conda_create("mofa_env")

#use it
use_condaenv("mofa_env", required = TRUE)

#install with pip
py_install("mofapy2", pip = TRUE)

#testing
reticulate::py_module_available("mofapy2")

#####################################

library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)
library(reticulate)

#loading data
utils::data("CLL_data")       
#checking dimensions ie number of features and samples
lapply(CLL_data,dim)
#readfile
CLL_metadata <- fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/cll_vignette/sample_metadata.txt")

#create mofa object
MOFAobject <- create_mofa(CLL_data)
#inspect ie summary table
MOFAobject

plot_data_overview(MOFAobject)

data_opts <- get_default_data_options(MOFAobject)
data_opts

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 15

model_opts

train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42

train_opts

MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts
)

#to run mofa, we need using mofapy
#conda create -n mofa_env python=3.8
#conda activate mofa_env
#pip install mofapy2

use_condaenv("mofa_env", conda = "C:/Users/asalehzade/AppData/Local/anaconda3/Scripts/conda.exe" ,required = TRUE)

MOFAobject <- run_mofa(MOFAobject)

#fixing non-numeric data issue
lapply(CLL_data, class)
lapply(CLL_data, function(x) class(x[1,1]))
lapply(CLL_data, colnames)

CLL_data <- lapply(CLL_data, function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x
})

#train then run
MOFAobject <- create_mofa(CLL_data)

MOFAobject <- prepare_mofa(MOFAobject)

MOFAobject <- run_mofa(MOFAobject, use_basilisk = TRUE)

########
slotNames(MOFAobject)

names(MOFAobject@data)

samples_metadata(MOFAobject) <- CLL_metadata

plot_factor_cor(MOFAobject)

plot_variance_explained(MOFAobject, max_r2=15)

plot_variance_explained(MOFAobject, plot_total = T)[[2]]

correlate_factors_with_covariates(MOFAobject, 
                                  covariates = c("Gender","died","age"), 
                                  plot="log_pval"
)

plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "Factor1"
)

plot_weights(MOFAobject,
             view = "Mutations",
             factor = 1,
             nfeatures = 5,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)

plot_top_weights(MOFAobject,
                 view = "Mutations",
                 factor = 1,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)

#distribution by frequency, wide is high
plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "IGHV",
            add_violin = TRUE,
            dodge = TRUE
)

plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "Gender",
            dodge = TRUE,
            add_violin = TRUE
)

plot_data_scatter(MOFAobject, 
                  view = "mRNA",
                  factor = 1,  
                  features = 4,
                  sign = "positive",
                  color_by = "IGHV"
) + labs(y="RNA expression")

plot_data_heatmap(MOFAobject, 
                  view = "mRNA",
                  factor = 1,  
                  features = 25,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)
