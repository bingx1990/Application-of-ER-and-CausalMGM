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
library(readr)

data <- read_delim("../dataset/age/age_all_days_imputed_5NN.txt", " ",
                   escape_double = FALSE, trim_ws = TRUE)


Y <- data$age
X <- data[,-(1:2)]

feature_names <- names(X)

################################################################################
#######                       ER on the whole dataset                    #######
################################################################################

X <- scale(X, T, T)
Y <- Y - mean(Y)
# 
# set.seed(20200316)
# delta_grid <- seq(0.35, 0.75, 0.01)
# 
# res <- ER(Y, X, delta_grid, pred = F, beta_est = "NULL", rep_CV = 50, verbose = T)

# res$K
# selected delta = 0.38

delta = 0.38
lbd = 0.5

res <- ER(Y, X, delta = delta, pred = F, beta_est = "LS", rep_CV = 50, verbose = T,
          merge = F, CI = T, correction = NULL, lbd = lbd)
res$K
### extract estimated beta, CIs and p-vals

# extract beta_hat and its standardized version
est_beta_LS <- res$beta
est_beta_LS_scaled <- diag(sqrt(diag(res$C))) %*% est_beta_LS
# return the ordered ranking based on the absolute values
order_coef_beta_LS <- order(abs(est_beta_LS_scaled), decreasing = T)
# 45  1 19 41  2 49 46 15  7 24 44 17 42  5 53  9 52 36

round(res$beta_CIs, 3)
p_vals <- 2 * pnorm(abs(res$beta / sqrt(res$beta_var / length(Y))), lower.tail = F)



res_dz <- ER(Y, X, delta = delta, pred = F, beta_est = "Dantzig", rep_CV = 50, verbose = T,
          merge = F, CI = T, correction = NULL, lbd = lbd)
round(res_dz$beta, 3)
Z_ind <- which(res_dz$beta != 0)

# extract beta_hat and its standardized version
est_beta_DZ <- res_dz$beta
est_beta_DZ_scaled <- diag(sqrt(diag(res_dz$C))) %*% est_beta_DZ
# return the ordered ranking based on the absolute values
order_coef_beta_DZ <- order(abs(est_beta_DZ_scaled), decreasing = T)


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
# 
# save_LOVE_result(list_clusters, delta, lbd, F, "age_dataset")
# save_LOVE_result(list_clusters_pure, delta, lbd, T, "age_dataset")


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


# file_name <- "../output/age_dataset/all_day_var_selection.txt"
# write(feature_names[coef_ind], file = file_name)
# write("\n", file = file_name, append = T)
# write(feature_names[Lasso_coef_ind], file = file_name, append = T)

######  2D Projection of X

Y <- data$age
Y_label <- rep(0, length(Y))
Y_label[Y >= 50] = 1

Z_hat <- X %*% A_hat %*% solve(crossprod(A_hat))
Z_tilde <- Pred_Z_BLP(X, A_hat, res$C, res$Gamma, 1:res$K)
Z_composite <- X[,coef_ind]
Z_Lasso <- X[,Lasso_coef_ind]

Plot_2D(Y_label, Z_hat, c("Elderly", "Adults"), axis_names = c("ER1", "ER2"))
Plot_2D(Y_label, Z_tilde[,Z_ind], c("Elderly", "Adults"), axis_names = c("ER1", "ER2"))
Plot_2D(Y_label, Z_composite, c("Elderly", "Adults"))
Plot_2D(Y_label, Z_Lasso, c("Elderly", "Adults"))


### Save the X matrix, predicted Z matrix and the estimated A

# 1  2  5  7  9 15 17 19 24 36 41 42 44 45 46 49 52 53
# write.csv(Z_tilde, file = "../output/age_dataset/Z_mat.csv", row.names = F)
# write.csv(A_hat, file = "../output/age_dataset/A_mat.csv", row.names = F)
# write.csv(cbind(Y, X), file = "../output/age_dataset/data_mat.csv", row.names = F)
# write.csv(X[,coef_ind], file = "../output/age_dataset/X_CR.csv", row.names = F)


################################################################################
######                       Prediction performance                       ######
################################################################################
X <- scale(X, T, T)

set.seed(20200309)

methods <- c("ER", "ER-Lasso-Dantzig", "PFR", "PLS", "Lasso", "Perm-Lasso")
# methods <- c("PFR", "PLS", "Lasso", "Perm-Lasso")

# result <- replicate(50, KfoldCV(10, Y = Y, X = X, methods = methods,
#                                 pt_option = "perm_X_Y", metric = "mse", standardize = F,
#                                 delta_grid = 0.38, ER_sparse = F))

result <- replicate(50, KfoldCV(5, Y = Y, X = X, methods = methods,
                                pt_option = "perm_X_Y", metric = "mse",
                                delta_grid = seq(0.35, 0.45, 0.02), ER_sparse = F))

CV_result <- t(sapply(result['MSE',], function(v) v))
fitted_array <- sapply(result['fitted',], function(v) v, simplify = "array")

fitted_mat <- apply(fitted_array, 1:2, mean)
CV_result <- CV_result / mean(Y ** 2)


# write.table(CV_result, file = "../output/age_dataset/all_days_CV_result_03_16_mse_sparse.txt", row.names = F, col.names = T)
# write.table(fitted_mat, file = "../output/age_dataset/all_days_CV_result_03_16_fit_sparse.txt", row.names = F, col.names = T)

# write.table(CV_result, file = "../output/age_dataset/all_days_CV_result_04_11_mse.txt", row.names = F, col.names = T)
# write.table(fitted_mat, file = "../output/age_dataset/all_days_CV_result_04_11_fit.txt", row.names = F, col.names = T)


# CV_result <- read.table("../output/age_dataset/all_days_CV_result_03_16_mse_sparse.txt", header = T)
# fitted_mat <- read.table("../output/age_dataset/all_days_CV_result_03_16_fit_sparse.txt", header = T)


CV_result <- read.table("../output/age_dataset/all_days_CV_result_04_11_mse.txt", header = T)
fitted_mat <- read.table("../output/age_dataset/all_days_CV_result_04_11_fit.txt", header = T)


# wilcox.test(CV_result$ER, CV_result$Perm.Lasso, alternative = "less")  ##  p-value < 2.2e-16

# sub_CV_result <- data.frame(sapply(1:ncol(CV_result), FUN = function(k, data) {
#   col_k <- data[,k]
#   col_k[which(col_k <= quantile(col_k, probs = 0.95))]
# }, data = CV_result))
# 
# names(sub_CV_result) <- names(CV_result)

methods[2] <- "Composite"

Violin_plot(methods, 1 - CV_result, MSE_flag = F, ylab_name = "1 - (normalized) mean square error",
            label_names = methods[c(1,2,5,3,4,6)], ylims = c(0.6, 1))

Fitted_plot(Y, fitted_mat, methods, order = c(1,2,5,3,4,6))


