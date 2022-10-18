
rm(list=ls())
setwd("../Code for ER")
source("SupLOVE.R")
setwd("..")
source("K-CV.R")
library(readr)
library(ROCR)
library(e1071)
library(doParallel)
library(foreach)
source("Helper.R")
source("Other_algorithms.R")



ma_data <- read_csv("total_data.csv")


raw_data <- ma_data[,-c(1,2)]

sub_data <- subset(raw_data, raw_data$time_point %in% c("T0", "T4"))
sub_data$time_point <- factor(sub_data$time_point, levels = c("T0", "T4"),
                              labels = c(0, 1))

Y <- as.numeric(levels(sub_data$time_point))[sub_data$time_point]
X <- data.frame(sub_data[,-1])

# levels(Y) <- c(0, 1, 2, 3, 4)

################################################################################
#######                       ER on the whole dataset                    #######
################################################################################

# X <- scale(X, T, T)
# Y <- Y - mean(Y)
# 
# set.seed(20200316)
# delta_grid <- seq(0.1, 0.2, 0.01)
# 
# res <- ER(Y, X, delta_grid, pred = F, beta_est = "NULL", rep_CV = 50, verbose = T, merge = F, parallel = T)
# 
# res$K   # 674
# # # selected delta = 0.14

# lbd = 0.5
# delta = 0.14
# 
# res <- ER(Y, X, delta = delta, pred = F, beta_est = "Dantzig", rep_CV = 50, verbose = T,
#           merge = F, CI = T, correction = "Bonferroni", lbd = lbd)
# res$K
# length(which(res$beta!=0))   # 86

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
# 



##### K-fold


X <- scale(X, T, T)

set.seed(20200320)

methods = c("ER", "ER-Lasso-Dantzig", "PLS", "PFR", "Lasso", "Perm-Lasso")

delta_grid <- 0.14

no_cores <- 10
registerDoParallel(cores=no_cores) 
cl <- makeCluster(no_cores)

result <- foreach(i = 1:50, .combine = "rbind") %dopar%
  KfoldCV_tenary(10, Y, X, methods, Y_type = "factor", pt_option = "perm_Y_before_split", 
                 delta_grid = delta_grid, ER_sparse = F)
stopCluster(cl)  

CV_result = t(sapply(result[,'metric'], function(v) v))
# pred_array <- sapply(result[,'pred'], function(v) v, simplify = "array")

# write.table(CV_result, file = "CV_result_T0T4.txt", row.names = F, col.names = T)


CV_result = read.table("CV_result_T0T4.txt", header = T)

Violin_plot(methods, CV_result, MSE_flag = F, label_names = methods[c(1,2,5,3,4,6)], ylims = c(0.2, 0.9))




# 
# 
# set.seed(20191205)
# 
# methods = c("PFR", "PLS", "Lasso", "Perm-Lasso")
# 
# result <- replicate(1, KfoldCV_tenary(60, Y, X, methods, Y_type = "count",
#                                       pt_option = "perm_Y_before_split"))
# CV_result <- t(sapply(result['metric',], function(v) v))
# pred_array <- sapply(result['pred',], function(v) v, simplify = "array")[,,1]
# 
# 
# true_groups <- list(which(Label_Y(Y) == 0), which(Label_Y(Y) == 1), which(Label_Y(Y) == 2))
# 
# cf_list <- c()
# 
# for (i in 1:length(methods)) {
#   pred_i <- pred_array[,i]
#   est_groups_i <- list(which(pred_i == 0), which(pred_i == 1), which(pred_i == 2))
#   cf_list[[i]] <- ConfusionMat(est_groups_i, true_groups)
# }
# 
# names(cf_list) <- methods
# 
# 
# ConfusionMat <- function(est_groups, true_groups = list(1:20, 21:40, 41:60)) {
#   col_count <- c()
#   for (i in 1:length(est_groups)) {
#     groupi <- est_groups[[i]]
#     col_count_i <- rep(0, 3)
#     for (j in groupi) {
#       true_ind <- which(sapply(1:3, FUN = function(x) {
#         j %in% true_groups[[x]]
#       }))
#       col_count_i[true_ind] = col_count_i[true_ind] + 1
#     }
#     col_count <- rbind(col_count, col_count_i)
#   }
#   rownames(col_count) <- c("Est G1", "Est G2", "Est G3")
#   colnames(col_count) <- c("True G1", "True G2", "True G3")
#   return(col_count)
# }

