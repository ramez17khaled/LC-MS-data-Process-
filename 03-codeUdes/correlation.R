## Heat map and discriminants features correlation between H and T
#################################################################################################################################################################################
dataPOS <- read.csv("D:/data/MSDial/04-codeOutput/Thermo_results/POS-manipulated.csv", sep = ";")
dataNEG <- read.csv("D:/data/MSDial/04-codeOutput/Thermo_results/NEG-manipulated.csv", sep=";")
volcano_H_T_dataPOS <- read.csv("D:/data/MSDial/04-codeOutput/Thermo_results/volcano_T_H_POS.csv")
volcano_H_T_dataNEG <- read.csv("D:/data/MSDial/04-codeOutput/Thermo_results/volcano_T_H_NEG.csv")

significant_metabolites_POS <- volcano_H_T_dataPOS[volcano_H_T_dataPOS$P_Value < 0.05, ]
significant_metabolites_NEG <- volcano_H_T_dataNEG[volcano_H_T_dataNEG$P_Value < 0.05, ]

# Find common metabolite names in data
M_columns_POS <- grep("^M", colnames(dataPOS), value = TRUE)
M_columns_NEG <- grep("^M", colnames(dataNEG), value = TRUE)

common_metabolites_POS <- intersect(names(dataPOS)[names(dataPOS) %in% M_columns_POS], significant_metabolites_POS$Metabolite)
common_metabolites_NEG <- intersect(names(dataNEG)[names(dataNEG) %in% M_columns_NEG], significant_metabolites_NEG$Metabolite)

# Filter columns in data based on common metabolite names
fdata_POS <- dataPOS[, c(common_metabolites_POS, "sample_name", "SampleType", "Condition", "Etude")]
fdata_NEG <- dataNEG[, c(common_metabolites_NEG, "sample_name", "SampleType", "Condition", "Etude")]

# num matrix
X_POS <- fdata_POS[, !(names(fdata_POS) %in% c("Condition", "SampleType", "sample_name", "Etude"))]
X_NEG <- fdata_NEG[, !(names(fdata_NEG) %in% c("Condition", "SampleType", "sample_name", "Etude"))]

# corr heatmap POS
corr_mat_POS <- round(cor(X_POS),2)  
head(corr_mat_POS)
library(reshape2)
melted_corr_mat_POS <- melt(corr_mat_POS)
library(heatmaply)
heatmaply_cor(x = cor(X_POS), xlab = "Features", 
              ylab = "Features", k_col = 2, k_row = 2)

# corr heatmap NEG
corr_mat_NEG <- round(cor(X_NEG),2)  
head(corr_mat_NEG)
library(reshape2)
melted_corr_mat_NEG <- melt(corr_mat_NEG)
library(heatmaply)
heatmaply_cor(x = cor(X_NEG), xlab = "Features", 
              ylab = "Features", k_col = 2, k_row = 2)

# Save correlation matrices as CSV files
write.csv(corr_mat_POS, "D:/data/MSDial/04-codeOutput/Thermo_results/correlation_matrix_H_T_POS.csv", row.names = TRUE)
write.csv(corr_mat_NEG, "D:/data/MSDial/04-codeOutput/Thermo_results/correlation_matrix_H_T_NEG.csv", row.names = TRUE)


