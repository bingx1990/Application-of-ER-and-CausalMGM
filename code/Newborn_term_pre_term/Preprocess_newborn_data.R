### Data preprocessing
rm(list=ls())
setwd("../../dataset/newborn")

library(readr)
Cell <- data.frame(read_csv("Cell population frequencies.csv"))
ComBat <- data.frame(read_csv("Short final_ComBat.csv"))

# View(Cell)
# View(ComBat)

# names(Cell)
# names(ComBat)


Cell <- Cell[,-c(5:6, 9:14)]
ComBat <- ComBat[,-c(1,3:15)]

All_data <- merge(Cell, ComBat, by = "Sample", all.x = T)
All_data$Time.point <- as.factor(All_data$Time.point)
All_data$Group <- as.factor(All_data$Group)
All_data$Relation <- as.factor(All_data$Relation)

dim(All_data)
# View(All_data)

All_data <- cbind(All_data, contain_NA = rowSums(is.na(All_data)) > 0)

# exclude the control samples ("Adult" in rows 326:337)
All_data <- All_data[1:325,]

# subset(All_data[which(All_data$contain_NA == T),], select = c("Subject", "Time.point", "Group", "contain_NA"))

All_data_no_father <- subset(All_data, All_data$Relation != "Father")
dim(All_data_no_father)   # 280 282


### For child

All_data_child_only <- subset(All_data, All_data$Relation == "Child")
dim(All_data_child_only)   # 183 282

child_NAs <- c("Cord blood" = mean(subset(All_data_child_only, All_data_child_only$Time.point == "Cord blood")$contain_NA),
  "Week 1" = mean(subset(All_data_child_only, All_data_child_only$Time.point == "Week 1")$contain_NA),
  "Week 4" = mean(subset(All_data_child_only, All_data_child_only$Time.point == "Week 4")$contain_NA),
  "Week 12" = mean(subset(All_data_child_only, All_data_child_only$Time.point == "Week 12")$contain_NA))

print(child_NAs, digits = 3)

### For mother

All_data_mother_only <- subset(All_data, All_data$Relation == "Mother")
dim(All_data_mother_only)  # 97 282

mother_NAs <- c("Birth" = mean(subset(All_data_mother_only, All_data_mother_only$Time.point == "Birth")$contain_NA),
                "Week 12" = mean(subset(All_data_mother_only, All_data_mother_only$Time.point == "Week 12")$contain_NA))

print(mother_NAs, digits = 3)



# save the data for "child"
# write.csv(All_data_child_only, "data_child.csv", row.names = F)


