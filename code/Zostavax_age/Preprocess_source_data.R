### Data preprocessing

rm(list=ls())
setwd("../../dataset/age/sourcedata")

library(readr)
age_data <- read_delim("VZV_age_gender2.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# View(age_data)

dim(age_data)  # 72 subjects 


### load feature datasets
load_data <- function(file_dir) {
  raw_data <- read_delim(file_dir, "\t", escape_double = FALSE, trim_ws = TRUE)
  new_data <- data.frame(t(raw_data[,-1]))
  colnames(new_data) <- raw_data[,1][[1]]
  new_data
}

data_list <- c()
count <- 1
for (file_i in dir()) {
  if (file_i != "ave_log2_QC_zerofil_vax05_neg_aplcms6_0828.txt" &
      file_i != "ave_log2_QC_zerofil_vax05_pos_aplcms6_0904.txt" &
      file_i != "VZV_age_gender2.txt") {
    data_list[[count]] <- load_data(file_i)
    names(data_list)[[count]] <- file_i
    count <- count + 1
    cat(file_i)
  }
}

# data_list$Tcell_IFNg3.txt

names(data_list)

dims <- c()
for (data_i in data_list) {
  dims <- rbind(dims, dim(data_i))
}








# ##### Extract data at Day 3
# # gene data at Day 3
# gene_data_D3 <- data_list[[8]][145:216,]
# rownames(gene_data_D3) <- sapply(rownames(gene_data_D3), function(x) strsplit(x,"_")[[1]][1], USE.NAMES = F)
# data_list[[8]] <- gene_data_D3
# 
# # BTMs data at Day 3
# BTMs_data_D3 <- data_list[[3]][,347:(2*346)]
# data_list[[3]] <- BTMs_data_D3
# 
# # MetaboClustersNeg and Pos at Day 3
# MetaNeg <- data_list[[4]][,39:76]
# data_list[[4]] <- MetaNeg
# 
# MetaPos <- data_list[[5]][,19:36]
# data_list[[5]] <- MetaPos
# 
# # Flow dataset at Day 3
# flow_data <- data_list[[9]]
# col_ind <- sapply(names(flow_data), function(x) {
#     ind1 <- strsplit(x,"_")[[1]][1] == "D3"
#     ind2 <- strsplit(x,"/")[[1]][1] == "D3"
#     return(ind1 | ind2)
#   }, USE.NAMES = F)
# data_list[[9]] <- flow_data[,col_ind]

# ##### Extract data at Day 7
# # gene data at Day 7
# gene_data_D7 <- data_list[[8]][217:288,]
# rownames(gene_data_D7) <- sapply(rownames(gene_data_D7), function(x) strsplit(x,"_")[[1]][1], USE.NAMES = F)
# data_list[[8]] <- gene_data_D7
# 
# # BTMs data at Day 7
# BTMs_data_D7 <- data_list[[3]][,(2*346+1):1038]
# data_list[[3]] <- BTMs_data_D7
# 
# # MetaboClustersNeg and Pos at Day 7
# MetaNeg <- data_list[[4]][,77:114]
# data_list[[4]] <- MetaNeg
# 
# MetaPos <- data_list[[5]][,37:54]
# data_list[[5]] <- MetaPos
# 
# # Flow dataset at Day 7
# flow_data <- data_list[[9]]
# col_ind <- sapply(names(flow_data), function(x) {
#   ind1 <- strsplit(x,"_")[[1]][1] == "D7"
#   ind2 <- strsplit(x,"/")[[1]][1] == "D7"
#   return(ind1 | ind2)
# }, USE.NAMES = F)
# data_list[[9]] <- flow_data[,col_ind]


# data indices that can be used 
data_ind <- which(dims[,1] > 36)


# exclude the gene dataset (the 8th)
data_ind <- setdiff(data_ind, 8)

### For data using only D3, we also exclude the BTM D0 dataset (the 2th)
data_ind <- setdiff(data_ind, 2)




names(data_list)[data_ind]

used_data_list <- lapply(data_list[data_ind], function(x) x)


####  unique subject id (with response)
ids <- sapply(age_data$ZV, function(x) {paste(strsplit(x, "-")[[1]], collapse = "")}, USE.NAMES = F)


#  unique subject id (with features)
ids_feature <- c()
for (i in 1:length(used_data_list))
  ids_feature <- union(ids_feature, rownames(used_data_list[[i]]))


# mismatched subject ids
ids_feature[which(ids_feature %in% ids == 0)]
ids[which(ids %in% ids_feature == 0)]

final_sub_ids <- intersect(ids, ids_feature)


# check how many missing ids for each dataset 

miss_ids <- list(length = length(used_data_list))
numb_nas <- c()
dims <- c()
for (i in 1:length(used_data_list)) {
  data_i <- used_data_list[[i]]
  id_ind_i <- which(rownames(data_i) %in% final_sub_ids)
  sub_data_i <- data_i[id_ind_i, ]
  used_data_list[[i]] <- sub_data_i
  
  dims <- rbind(dims, dim(sub_data_i))
  miss_id_ind <- which(final_sub_ids %in% rownames(sub_data_i) == 0)
  if (length(miss_id_ind) == 0) 
    miss_ids[[i]] <- 0
  else 
    miss_ids[[i]] <- miss_id_ind
  numb_nas[i] <- sum(is.na(sub_data_i))
}
rownames(dims) <- names(miss_ids) <- names(numb_nas) <- names(used_data_list)

dims
miss_ids


sub_age_data <- subset(age_data,ids %in% final_sub_ids)
data_mat <- data.frame(id = ids[ids %in% final_sub_ids], age = sub_age_data$AGE,
                       row.names = ids[ids %in% final_sub_ids])

for (i in 1:length(used_data_list)) {
  data_i <- used_data_list[[i]] 
  data_mat <- merge(data_mat, data.frame(id = rownames(data_i), data_i), by.x = "id",
                        all.x = T)
}

dim(data_mat)
# write.table(data_mat, "../raw_age_data_all_days.txt", row.names = T, col.names = T)




##### Inpute NAs by 5-NN 

raw_data <- read.table("../raw_age_data_all_days.txt")


raw_X <- raw_data[,-(1:2)]
data_matrix <- t(raw_X)
dim(data_matrix)    # 431 * 67
sum(is.na(data_matrix))   # 245 missing values
na_col_ind <- which(is.na(colSums(data_matrix)))
# each subject contains NAs miss 7:32 rows
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

imputed_data <- Impute_K_NN(5, na_col_ind, data_matrix)

X = t(imputed_data)
whole_data <- data.frame(id = raw_data$id, age = raw_data$age, X)
# write.table(whole_data,
#             file = "../age_all_days_imputed_5NN.txt",
#             col.names = T, row.names = F)
