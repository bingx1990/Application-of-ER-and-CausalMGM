
rm(list=ls())
setwd("../Code for ER")
source("SupLOVE.R")
setwd("../")
source("K-CV.R")
library(readr)
library(ROCR)
library(e1071)
library(doParallel)
library(foreach)
source("Helper.R")
source("Other_algorithms.R")
setwd("./malaria-code/")


# #### Load and preprocess the dataset
# library(readr)
# GPL571 <- read_delim("../dataset/malaria_data/GSE18323-GPL571_series_matrix.txt", 
#                      "\t", escape_double = FALSE, trim_ws = TRUE, 
#                      skip = 68)
# GPL570 <- read_delim("../dataset/malaria_data/GSE18323-GPL570_series_matrix.txt", 
#                      "\t", escape_double = FALSE, trim_ws = TRUE, 
#                      skip = 69)
# 
# 
# # View(GPL571)
# GPL571 <- GPL571[-nrow(GPL571),]
# GPL570 <- GPL570[-nrow(GPL570),]
# 
# X_571 <- data.frame(t(GPL571[,-1]))
# names(X_571) <- GPL571[,1][[1]]
# 
# X_570 <- data.frame(t(GPL570[,-1]))
# names(X_570) <- GPL570[,1][[1]]
# 
# # use the common gene-expressions of both dataset
# 
# X_570_sub <- X_570[,names(X_570) %in% names(X_571)]
# 
# dim(X_570_sub)
# dim(X_571)
# 
# 
# # the design matrix
# X_total <- rbind(X_570_sub, X_571)
# X_total <- cbind(test_id = rownames(X_total), X_total)
# 
# response <- read_delim("../dataset/malaria_data/response.csv", 
#                        "\t", escape_double = FALSE, trim_ws = TRUE, col_names = F)
# time_mat <- t(sapply(response$X2, FUN = function(x) {
#   unlist(strsplit(x, "    "))
# }))
# 
# id_mat <- data.frame(test_id = response$X1, subject_id = time_mat[,1], time_point = time_mat[,2])
# rownames(id_mat) <- NULL
# 
# 
# sub_X_total <- subset(X_total, X_total$test_id %in% id_mat$test_id)
# 
# total_data <- merge(id_mat, sub_X_total)
# 
# ## save the data=
# write_csv(total_data, "../dataset/malaria_data/total_data.csv")



ma_data <- read_csv("total_data.csv")


raw_data <- ma_data[,-c(1,2)]

sub_data <- subset(raw_data, raw_data$time_point %in% c("T0", "T2", "T4"))
sub_data$time_point <- factor(sub_data$time_point, levels = c("T0", "T2", "T4"),
                              labels = c(0, 1, 2))

Y <- as.numeric(levels(sub_data$time_point))[sub_data$time_point]
X <- data.frame(sub_data[,-1])

# levels(Y) <- c(0, 1, 2, 3, 4)

################################################################################
#######                       ER on the whole dataset                    #######
################################################################################

X <- scale(X, T, T)
Y <- Y - mean(Y)

# set.seed(20200316)
# delta_grid <- seq(0.01, 0.15, 0.01)
# 
# res <- ER(Y, X, delta_grid, pred = F, beta_est = "NULL", rep_CV = 50, verbose = T, merge = F, parallel = T)
# 
# res$K
# selected delta = 0.04

lbd = 0.5
delta = 0.04

res <- ER(Y, X, delta = delta, pred = F, beta_est = "Dantzig", rep_CV = 50, verbose = T,
          merge = F, CI = T, correction = "Bonferroni", lbd = lbd)
res$K
length(which(res$beta != 0))   # 86

# save.image("ER_result.RData")

load("ER_result.RData")

# clusters <- recoverGroup(res$A)
# 
# sapply(clusters, function(x) as.vector(sort(unlist(x))))
# sapply(res$I_ind, function(x) as.vector(unlist(x)))

#
# # output cluster results
# # save_LOVE_result <- function(clusters, delta, lbd, only.pure = F) {
# #   K <- length(clusters)
# #   if (only.pure)
# #     savePath = paste(format(Sys.time(),"../output/Cibersort_dataset/love_%m.%d_nGroup"),K,"delta=",delta,"lbd=",lbd,"_pure",sep='')
# #   else
# #     savePath = paste(format(Sys.time(),"../output/Cibersort_dataset/love_%m.%d_nGroup"),K,"delta=",delta,"lbd=",lbd, " ", sep='')
# #   dir.create(file.path(getwd(),savePath))
# #
# #   foreach (k = 1:K) %do% {
# #     groupk = sort(unlist(clusters[[k]], use.names = F))
# #     name = file.path(getwd(), savePath, paste("group",k,"-size-",length(groupk),
# #                                             ".txt", sep=''))
# #     write.table(groupk, file=name, sep="\t", col.name=FALSE, row.names=FALSE)
# #   }
# #   return(savePath)
# # }
#
# save_LOVE_result(clusters, delta, lbd)
# save_LOVE_result(res$I_ind, delta, lbd, T)
#


Z_ind <- which(res$beta != 0)

A_hat <- res$A
feature_ind <- which(rowSums(abs(A_hat[,Z_ind, drop = F])) != 0)

set.seed(20200415)

cvfit = cv.glmnet(X[,feature_ind,drop = F], Y, standardize = F, intercept = F)
sign_coef <- predict(cvfit, s = "lambda.min", type = "coefficient")[-1]
coef_ind <- feature_ind[which(sign_coef != 0)]


Lasso_cvfit = cv.glmnet(X, Y, standardize = F, intercept = F)
lbd_ind <- which(Lasso_cvfit$lambda == Lasso_cvfit$lambda.min)
Lasso_sign_coef <- predict(Lasso_cvfit, s = Lasso_cvfit$lambda[lbd_ind], type = "coefficient")[-1]
Lasso_coef_ind <- which(Lasso_sign_coef != 0)

cat("ER-Lasso selects", length(coef_ind), "features while Lasso selects",
    length(Lasso_coef_ind), "features with", length(intersect(Lasso_coef_ind, coef_ind)),
    "features in common.")

# feature_names <- colnames(X)
# file_name <- "../../output/malaria_dataset/var_selection.txt"
# write(feature_names[coef_ind], file = file_name)
# write("\n", file = file_name, append = T)
# write(feature_names[Lasso_coef_ind], file = file_name, append = T)


# ######  2D Projection of X

Y_label <- as.numeric(levels(sub_data$time_point))[sub_data$time_point]

Z_hat <- X %*% A_hat %*% solve(crossprod(A_hat))
Z_tilde <- Pred_Z_BLP_avg(X, A_hat, res$C, 1:res$K)
Z_composite <- X[,coef_ind]
Z_Lasso <- X[,Lasso_coef_ind]

Plot_2D(Y_label, Z_hat[,Z_ind], c("Pre", "Post1", "Post2"), axis_names = c("ER1", "ER2"))
Plot_2D(Y_label, Z_composite, c("Pre", "Post1", "Post2"))
Plot_2D(Y_label, Z_Lasso, c("Pre", "Post1", "Post2"))


sub_Y_label <- Y_label %in% c(0, 1)
Plot_2D(Y_label[sub_Y_label], Z_hat[sub_Y_label,Z_ind], c("Pre", "Post1"), axis_names = c("ER1", "ER2"))
Plot_2D(Y_label[sub_Y_label], Z_composite[sub_Y_label,], c("Pre", "Post1"))

sub_Y_label <- Y_label %in% c(0, 2)
new_Y_label <- rep(0, length(Y_label))
new_Y_label[Y_label == 2] = 1
Plot_2D(new_Y_label[sub_Y_label], Z_hat[sub_Y_label,Z_ind], c("Pre", "Post2"), axis_names = c("ER1", "ER2"))
Plot_2D(new_Y_label[sub_Y_label], Z_composite[sub_Y_label,], c("Pre", "Post2"))

sub_Y_label <- Y_label %in% c(1, 2)
Plot_2D(Y_label[sub_Y_label], Z_hat[sub_Y_label,Z_ind], c("Post1", "Post2"), axis_names = c("ER1", "ER2"))
Plot_2D(Y_label[sub_Y_label], Z_composite[sub_Y_label,], c("Post1", "Post2"))



### Save the X matrix, predicted Z matrix and the estimated A

### Z_ind:  36   64   85  165  185  206  211  212  221  230  251  263  267  270  297  302  314
#  319  323  331  333  342  344  369  422  459  481  514  558  561  568  574  580  641
#  661  703  765  786  812  815  847  898  928  985 1006 1025 1048 1071 1084 1115 1128
# 1166 1176 1179 1200 1210 1262 1294 1327 1340 1349 1361 1370 1377 1383 1428 1434 1439
# 1442 1481 1483 1493 1500 1513 1515 1534 1539 1542 1559 1566 1573 1590 1614 1641 1646
# 1653
write.csv(Z_tilde, file = "../output/malaria_dataset/Z_mat.csv", row.names = F)
write.csv(A_hat, file = "../output/malaria_dataset/A_mat.csv", row.names = F)
write.csv(cbind(Y, X), file = "../output/malaria_dataset/data_mat.csv", row.names = F)

write.csv(X[,coef_ind], file = "../output/malaria_dataset/X_CR.csv", row.names = F)

############################################################
#####                      K-fold                     ######
############################################################


X <- scale(X, T, T)

set.seed(20200320)

methods = c("ER", "ER-Lasso-Dantzig", "PLS", "PFR", "Lasso", "Perm-Lasso")
# methods = "ER-Lasso-Dantzig"

delta_grid <- 0.04

no_cores <- 10
registerDoParallel(cores=no_cores) 
cl <- makeCluster(no_cores)

result <- foreach(i = 1:50, .combine = "rbind") %dopar%
  KfoldCV_tenary(10, Y, X, methods, Y_type = "factor", pt_option = "perm_Y_before_split", 
                 delta_grid = delta_grid, ER_sparse = F)
stopCluster(cl)  

CV_result = t(sapply(result[,'metric'], function(v) v))
# pred_array <- sapply(result[,'pred'], function(v) v, simplify = "array")

# write.table(CV_result, file = "CV_result.txt", row.names = F, col.names = T)


methods = c("ER", "Composite", "PLS", "PFR", "Lasso", "Perm-Lasso")

CV_result = read.table("CV_result.txt", header = T)

Violin_plot(methods, CV_result, MSE_flag = F, label_names = methods[c(1,2,5,3,4,6)], ylims = c(0.2, 0.9))




#### LOOCV 

X <- scale(X, T, T)

set.seed(20191205)

methods = c("ER", "ER-Lasso-Dantzig", "PFR", "PLS", "Lasso", "Perm-Lasso")

# result <- KfoldCV_tenary(10, Y, X, methods, Y_type = "factor", pt_option = "perm_Y_before_split", delta_grid = 0.04)
# pred_array <- result$pred


pred_array <- KfoldCV_tenary_parallel(nrow(X), Y, X, methods, Y_type = "factor", pt_option = "perm_Y_before_split", delta_grid = 0.04)

# write.table(pred_array, file = "pred_result_CV10.txt", row.names = F, col.names = F)
# write.table(pred_array, file = "pred_result_LOOCV.txt", row.names = F, col.names = F)

pred_array <- read.table("pred_result_LOOCV.txt")

ConfusionMat <- function(est_groups, true_groups = list(1:20, 21:40, 41:60)) {
  col_count <- c()
  for (i in 1:length(est_groups)) {
    groupi <- est_groups[[i]]
    col_count_i <- rep(0, 3)
    for (j in groupi) {
      true_ind <- which(sapply(1:3, FUN = function(x) {
        j %in% true_groups[[x]]
      }))
      col_count_i[true_ind] = col_count_i[true_ind] + 1
    }
    col_count <- rbind(col_count, col_count_i)
  }
  rownames(col_count) <- c("Est G1", "Est G2", "Est G3")
  colnames(col_count) <- c("True G1", "True G2", "True G3")
  return(col_count)
}

true_groups <- list(which(Y == 0), which(Y == 1), which(Y == 2))

cf_list <- c()

for (i in 1:length(methods)) {
  pred_i <- pred_array[,i]
  est_groups_i <- list(which(pred_i == 1), which(pred_i == 2), which(pred_i == 3))
  cf_list[[i]] <- ConfusionMat(est_groups_i, true_groups)
}

names(cf_list) <- methods

sapply(cf_list, function(x) sum(diag(x)) / length(Y))





library(gridExtra)
library(grid)
cm_ER <- cf_list$ER
n <- length(Y)

tt <- ttheme_minimal(base_size = 12,
                     core=list(bg_params = list(fill = blues9[(rank(cm_ER) + 3) / 2], col=NA),
                               fg_params=list(fontface=1), alpha = 1),
                     colhead=list(fg_params=list(col="blue", fontface= 1)),
                     rowhead=list(fg_params=list(col="blue", fontface= 1)))

g_ER <- tableGrob(round(cm_ER / n, 3), theme = tt)
g_CR <- tableGrob(round(cf_list$Composite / n, 3), theme = tt)
g_Lasso <- tableGrob(round(cf_list$Lasso / n, 3), theme = tt)
g_PLS <- tableGrob(round(cf_list$PLS / n, 3), theme = tt)
g_PFR <- tableGrob(round(cf_list$PFR / n, 3), theme = tt)
grid.newpage()
grid.arrange(arrangeGrob(g_ER, top = "ER"), 
             arrangeGrob(g_CR, top = "Composite"), 
             arrangeGrob(g_Lasso, top = "Lasso"), 
             arrangeGrob(g_PLS, top = "PLS"), 
             arrangeGrob(g_PFR, top = "PFR"), nrow = 3)
