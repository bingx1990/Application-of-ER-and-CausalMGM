
rm(list=ls())
setwd("../Code for ER")
source("SupLOVE.R")
setwd("..")
source("K-CV.R")
source("Helper.R")
source("Other_algorithms.R")
library(readr)
library(ROCR)

child_data <- read.csv("../dataset/newborn/data_child.csv")
# View(child_data)

# ### Week 1    dims: 56  282 
# subdata <- subset(child_data, child_data$Time.point == "Week 1")

# ### Week 4    dims: 29 275
# subdata <- subset(child_data, child_data$Time.point == "Week 4")
# 
### Week 12   dims: 46 282
subdata <- subset(child_data, child_data$Time.point == "Week 12")


raw_Y <- subdata$Group
Y <- rep(0, length(raw_Y))
Y[raw_Y == "Premature"] = 1
raw_X <- subdata[,-c(1:7, 282)]   # the 282th column is the NA flag

feature_names <- names(subdata)[-c(1:7, 282)]

### Impute data via KNN
Impute_K_NN <- function(K, na_col_ind, data_matrix) {
  # na_row_ind <- 7:32
  new_data <- data_matrix
  for (sub_ind in na_col_ind) {
    na_row_ind <- which(is.na(data_matrix[, sub_ind]))
    sub_data <- data_matrix[-na_row_ind, sub_ind]
    K_neighbor <- CompKNN(K, sub_data, na_col_ind, na_row_ind, data_matrix)
    new_data[na_row_ind, sub_ind] <- apply(data_matrix[na_row_ind, K_neighbor], 1, mean)
  }
  return(new_data)
}
CompKNN <- function(K, sub_data, na_col_ind, na_row_ind, data_matrix) {
  candidates <- setdiff(1:ncol(data_matrix), na_col_ind)
  distance <- rep(Inf, length(candidates))
  distance[candidates] <- sapply(candidates, FUN = function(x) {
    mean((data_matrix[-na_row_ind, x] - sub_data) ** 2)
  })
  neighbor_ind <- order(distance)[1:K]
  return(neighbor_ind)
}

data_matrix <- t(raw_X)
na_col_ind <- which(is.na(colSums(data_matrix)))
imputed_data <- Impute_K_NN(5, na_col_ind, data_matrix)

X <- t(imputed_data)


################################################################################
#######                       ER on the whole dataset                    #######
################################################################################

X <- scale(X, T, T)
Y <- Y - mean(Y)

# set.seed(20200316)
# delta_grid <- seq(0.1, 0.35, 0.01)
#
# res <- ER(Y, X, delta_grid, pred = F, beta_est = "NULL", rep_CV = 100, verbose = T, merge = F)
#
# res$K
# # selected delta = 0.15

lbd = 1
delta = 0.15

res <- ER(Y, X, delta = delta, pred = F, beta_est = "LS", rep_CV = 50, verbose = T,
          merge = F, CI = T, correction = "Bonferroni", lbd = lbd)
res$K    # 14
round(res$beta_CIs, 3)
p_vals <- 2 * pnorm(abs(res$beta) / sqrt(res$beta_var / length(Y)), lower.tail = F)
round(p_vals, 3)
# 
# Z_ind <- which(p_vals <= 0.1)
# 
# clusters <- recoverGroup(res$A)
# 
# list_clusters_pure <- lapply(Z_ind, FUN = function(x, data_list, feature_name) {
#   feature_name[unlist(data_list[[x]])]
# }, data_list = res$I_ind, feature_name = colnames(X))
# 
# 
# list_clusters <- lapply(Z_ind, FUN = function(x, data_list, feature_name) {
#   feature_name[unlist(data_list[[x]])]
# }, data_list = clusters, feature_name = colnames(X))
# 
# 
# sapply(list_clusters_pure, FUN = function(x) {length(x)})
# sapply(list_clusters, FUN = function(x) {length(x)})
# 
# # save_LOVE_result(list_clusters, delta, lbd, F, "newborn_dataset/clusters_week12")
# # save_LOVE_result(list_clusters_pure, delta, lbd, T, "newborn_dataset/clusters_week12")




Z_ind <- which(p_vals <= 0.05)
A_hat <- res$A
feature_ind <- which(rowSums(abs(A_hat[,Z_ind, drop = F])) != 0)

set.seed(20200415)

cvfit = cv.glmnet(X[,feature_ind,drop = F], Y, standardize = F, intercept = F)
sign_coef <- predict(cvfit, s = "lambda.min", type = "coefficient")[-1]
coef_ind <- feature_ind[which(sign_coef != 0)]


Lasso_cvfit = cv.glmnet(X, Y, standardize = F, intercept = F)
Lasso_sign_coef <- predict(Lasso_cvfit, s = "lambda.min", type = "coefficient")[-1]
Lasso_coef_ind <- which(Lasso_sign_coef != 0)

cat("ER-Lasso selects", length(coef_ind), "features while Lasso selects",
    length(Lasso_coef_ind), "features with", length(intersect(Lasso_coef_ind, coef_ind)),
    "features in common.")


# file_name <- "../output/newborn_dataset/week12_var_selection.txt"
# write(feature_names[coef_ind], file = file_name)
# write("\n", file = file_name, append = T)
# write(feature_names[Lasso_coef_ind], file = file_name, append = T)


######  2D Projection of X

Y_label <- rep(0, length(raw_Y))
Y_label[raw_Y == "Premature"] = 1

Z_hat <- X %*% A_hat %*% solve(crossprod(A_hat))
Z_tilde <- Pred_Z_BLP(X, A_hat, res$C, res$Gamma, 1:res$K)
Z_composite <- X[,coef_ind]
Z_Lasso <- X[,Lasso_coef_ind]

Plot_2D(Y_label, Z_hat, c("Pre-term", "Term"))
Plot_2D(Y_label, Z_tilde[,Z_ind], c("Pre-term", "Term"), 
        axis_names = c("ER1", "ER2"))
Plot_2D(Y_label, Z_composite, c("Pre-term", "Term"), "topleft")
Plot_2D(Y_label, Z_Lasso, c("Pre-term", "Term"), "bottomleft")



### Save the X matrix, predicted Z matrix and the estimated A

# 5 7
write.csv(Z_tilde, file = "../output/newborn_dataset/Z_mat.csv", row.names = F)
write.csv(A_hat, file = "../output/newborn_dataset/A_mat.csv", row.names = F)
write.csv(cbind(Y, X), file = "../output/newborn_dataset/data_mat.csv", row.names = F)



################################################################################
######                    Classification performance                      ######
################################################################################
X <- scale(X, T, T)

methods = c("ER", "ER-Lasso-LS", "PFR", "PLS", "Lasso", "Perm-Lasso")

# delta <- seq(0.1, 0.3, 0.01)  # week 1
delta <- seq(0.1, 0.35, 0.01)   # week 12

set.seed(20200725)
pred_val <- CV_binary(Y, X, methods, delta, F)

write.table(pred_val, file = "../output/newborn_dataset/week12_non_sparse.txt", row.names = F, col.names = T)
# write.table(pred_val, file = "../output/newborn_dataset/week12_sparse.txt", row.names = F, col.names = T)

pred_val <- read.table("../output/newborn_dataset/week12_non_sparse.txt", header = T)
# pred_val <- read.table("../output/newborn_dataset/week12_sparse.txt", header = T)

methods[2] <- "Composite"

colnames(pred_val) <- methods

true_label <- Y
result <- Compute_metric(pred_val, true_label)

# auc
names(result$auc) <- methods
round(result$auc, 3)


### Week 12
## non-sparse 
# ER  Composite        PFR        PLS      Lasso     Perm-Lasso 
# 0.869      0.770      0.716      0.750      0.925      0.466 

## sparse
# ER         PFR        PLS        Lasso      Perm-Lasso 
# 0.817      0.716      0.750      0.921      0.601 

# Plot_ROC(result$error, methods, "bottomright")
Plot_ROC(result$error, setdiff(methods, "Perm-Lasso"), "bottomright")




################################################################################
######                 k-fold classification performance                  ######
################################################################################
library(e1071)
# standardize the features
X <- scale(X, T, T)

set.seed(20200320)

methods = c("ER", "ER-Lasso-LS", "PFR", "PLS", "Lasso", "Perm-Lasso")
# delta_grid <- seq(0.1, 0.3, 0.01)   # grid for week 1
delta_grid <- seq(0.1, 0.35, 0.01)   # grid for week 12

result <- replicate(50, KfoldCV_tenary(5, Y, X, methods, Y_type = "factor",
                                       pt_option = "perm_Y_before_split",
                                       delta_grid = delta_grid,
                                       ER_sparse = F))

CV_result <- t(sapply(result['metric',], function(v) v))
pred_array <- sapply(result['pred',], function(v) v, simplify = "array")
# summary(CV_result, digits = 3)

# write.table(CV_result, file = "../output/newborn_dataset/week12_accuracy_04_11.txt", row.names = F, col.names = T)
# write.table(CV_result, file = "../output/newborn_dataset/week12_accuracy_sparse.txt", row.names = F, col.names = T)

CV_result <- read.table("../output/newborn_dataset/week12_accuracy_04_11.txt", header = T)
# CV_result <- read.table("../output/newborn_dataset/week12_accuracy_sparse.txt", header = T)


methods[2] <- "Composite"

Violin_plot(methods, CV_result, MSE_flag = F, ylims = c(0.25, 1),
            label_names = methods[c(1,2,5,3,4,6)])











##################################################
######
######        check clusters overlap       #######
######
##################################################

rm(list=ls())
setwd("../output/newborn_dataset")
library(readr)

# # alpha = 0.1 level
# path_w1 <- "./clusters_week1/love_04.08_nGroup5delta=0.17lbd=1 /"
# path_w12 <- "./clusters_week12/love_04.08_nGroup3delta=0.15lbd=1 /"

# alpha = 0.05 level
path_w1 <- "./clusters_week1/love_04.04_nGroup5delta=0.17lbd=1/"
path_w12 <- "./clusters_week12/love_04.04_nGroup2delta=0.15lbd=1/"

clusters_w1 <- clusters_w12 <- c()

for (i in 1:length(dir(path_w1))) {
  clusters_w1[[i]] <- read_csv(paste(path_w1, dir(path_w1)[i], sep = ""), col_names = FALSE)
}  

for (i in 1:length(dir(path_w12))) {
  clusters_w12[[i]] <- read_csv(paste(path_w12, dir(path_w12)[i], sep = ""), col_names = FALSE)
} 


common <- matrix(0, length(clusters_w1), length(clusters_w12),
                 dimnames = list(sapply(clusters_w1, function(x) length(x$X1)),
                                 sapply(clusters_w12, function(x) length(x$X1))))

for (i in 1:length(clusters_w1)) {
  for (j in 1:length(clusters_w12)) {
    common[i,j] <- length(intersect(clusters_w1[[i]]$X1, clusters_w12[[j]]$X1))
  }
}

common

