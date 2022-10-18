########################################################################
#######                       Predictions                       ########
########################################################################
rm(list=ls())
setwd("Code for ER")
source("SupLOVE.R")
library(pls)

set.seed(20190317)
### Varying ps
n <- 300
ps <-  seq(200, 1000, 200)
K <- 10
m <- 5
Npure <- rep(list(c(m, 0)), each = K)
C = GenC(seq(2.5, 3, length.out = K), nu = 0.3)
beta <- runif(K, 0, 1)

rep_time <- 50
errors_NPR <- rank <- rank_pcr_ratio <- snr <- c()

for (p in ps) {
  errors_NPR_item <- c()
  rank_item <- rank_pcr_ratio_item <- c()
  A <- GenA(p, Npure, 10, 2, 0.5, "gaussian")
  # A <- GenA(p, Npure, maxNumb = K %/% 2 + 1, minNumb = 2)
  Gamma <- runif(p, 3, 5)
  snr <- c(snr, eigen(A %*% C %*% t(A), only.values = T)$values[K] / max(Gamma))
  
  for (i in 1:rep_time) {
    data <- Data_gen(n, p, beta, C, A, Gamma)
    X <- scale(data$X, T, F);  Y <- scale(data$Y, T, F)
    
    result_LOVE <- LOVE(X, delta = seq(0.1, 0.6, 0.1))
    I_hat <- result_LOVE$I
    A_hat <- result_LOVE$A
    rank_item[i] <- ncol(A_hat)
    
    result <- Pred_avg(Y, X, A_hat, C_hat, 1:ncol(X))
    
    result_PCR <- PCR(X, Y, K)
    
    result_PCR_ratio <- PCR(X, Y, option = "ratio")
    rank_pcr_ratio_item[i] <- result_PCR_ratio$K
    
    
    ### New data prediction
    new_data <- Data_gen(n, p, beta, C, A, Gamma)
    new_X <- new_data$X
    new_mean <- (new_data$Z) %*% beta 
    
    ### PLS 
    data = data.frame(Y = rbind(Y, new_mean), X = rbind(X, new_X))
    train_ind <- 1:n
    model_pls <- plsr(Y ~ ., ncomp = 15, data = data[train_ind, ], validation = "CV") 
    n_comp <- selectNcomp(model_pls, method = "randomization")
    n_comp <- ifelse(n_comp == 0, 1, n_comp)
    pred_pls <- predict(model_pls, comps = n_comp, newdata = data[-train_ind,], 
                        type = "response")
    
    NPRs <- c("Full-avg" = Mse(new_mean, new_X %*% result$theta),
              "PCR-oracle" = Mse(new_mean, new_X %*% result_PCR$theta),
              "PCR-ratio" = Mse(new_mean, new_X %*% result_PCR_ratio$theta),
              "Lasso" = Mse(new_mean, LASSO(X, Y, new_X)),
              "PLS" = Mse(new_mean, pred_pls))
    errors_NPR_item <- rbind(errors_NPR_item, NPRs)
    
    if (i %% 25 == 0 ) {
      cat("Finishing the i =", i, "th iteration for p = ",p,"...\n")
    }
  }
  rank <- c(rank, mean(rank_item))
  rank_pcr_ratio <- c(rank_pcr_ratio, mean(rank_pcr_ratio_item))
  errors_NPR <- rbind(errors_NPR, apply(errors_NPR_item, 2, mean))
}

result <- cbind(errors_NPR, rank, rank_pcr_ratio)

rownames(result) <- c("p = 200", "p = 400", "p = 600", "p = 800", "p = 1000")

# write.table(result, file = "../output/syn_mse.txt")

data_pred <- read.table("../output/syn_mse.txt")

# data_pred <- data_pred[,-4]
x_labs <- c(200, 400, 600, 800, 1000)



pred_names <- c("ER", "PCR", "Lasso", "PLS")

par(mgp = c(1.3,0.3,0), mar = c(2.5,1.3,1.5,1), cex.lab = 1, cex.axis = .9, tcl = -0.25)
plot(x_labs, data_pred[,1], col = 2, lty = 1, pch = 1, type = "b", xlab = "p", 
     ylab = "", 
     ylim = c(0, max(data_pred[,1:5])) + 0.5 * c(-1, 1), lwd = 1.5,
     main = "n = 300, K = 10, m = 5")
points(x_labs, data_pred[,3], col = 3, lty = 2, pch = 2, type = "b", lwd = 1.5)
points(x_labs, data_pred[,4], col = 4, lty = 3, pch = 3, type = "b", lwd = 1.5)
points(x_labs, data_pred[,5], col = 5, lty = 4, pch = 4, type = "b", lwd = 1.5)
legend("bottomleft", legend = pred_names, col = 2:5, lty = 1:4, 
       pch = 1:4, lwd = 1.5, seg.len = 3, cex = 1)




# pred_names <- c("ER", "PFR-oracle", "PFR")
# par(mgp = c(1.3,0.3,0), mar = c(2.5,1.3,1.5,1), cex.lab = 1, cex.axis = .9, tcl = -0.25)
# plot(x_labs, data_pred[,1], col = 2, lty = 1, pch = 1, type = "b", xlab = "p", 
#      ylab = "", 
#      ylim = range(data_pred[,1:length(pred_names)]) + 0.2 * c(-1, 1), lwd = 1.5,
#      main = "n = 300, K = 10, m = 5")
# points(x_labs, data_pred[,2], col = 4, lty = 3, pch = 3, type = "b", lwd = 1.5)
# points(x_labs, data_pred[,3], col = 5, lty = 4, pch = 4, type = "b", lwd = 1.5)
# legend("topright", legend = pred_names, col = c(2,4,5), lty = c(1,3,4), 
#        bty = "n", pch = c(1,3,4), lwd = 1.5, seg.len = 3, cex = 1)



