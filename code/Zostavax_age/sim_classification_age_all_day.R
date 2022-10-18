################################################################################
##########                                                            ##########
##########                      ER  application                      ###########
##########           prediction / classification  performance        ###########
##########                                                            ##########
################################################################################
rm(list=ls())
setwd("../Code for ER")
source("SupLOVE.R")
setwd("..")
source("K-CV.R")
source("Helper.R")
source("Other_algorithms.R")
library(ROCR)
library(readr)

data <- read_delim("../dataset/age/age_all_days_imputed_5NN.txt", " ",
                   escape_double = FALSE, trim_ws = TRUE)
Y <- data$age
X <- data[,-(1:2)]


################################################################################
######                    Classification performance                      ######
################################################################################
X <- scale(X, T, T)

methods = c("ER", "ER-Lasso-Dantzig", "PFR", "PLS", "Lasso", "Perm-Lasso")
methods = "PFR-3"

set.seed(20200316)

pred_val <- CV_binary(Y, X, methods, seq(0.35, 0.45, 0.02), F)

# write.table(pred_val, file = "../output/age_dataset/pred_val_binary_all_days.txt", row.names = F, col.names = T)

pred_val <- read.table("../output/age_dataset/pred_val_binary_all_days.txt", header = T)

methods[2] <- "Composite"

colnames(pred_val) <- methods

true_label <- rep(0, length(Y))
true_label[Y >= 50] = 1

result <- Compute_metric(pred_val, true_label)

# auc  
names(result$auc) <- methods
round(result$auc, 2)    

# ER  Composite        PFR        PLS      Lasso Perm-Lasso 
# 0.79       0.74       0.47       0.60       0.65       0.53  


Plot_ROC(result$error, methods)
Plot_ROC(result$error, setdiff(methods, "Perm-Lasso"), "bottomright")














