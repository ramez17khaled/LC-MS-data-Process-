

InputcsvFilePOS <- "D:/data/MSDial/04-codeOutput/Thermo_results/POS-manipulated.csv"
InputcsvFileNEG <- "D:/data/MSDial/04-codeOutput/Thermo_results/NEG-manipulated.csv"
OutputcsvFile <- "D:/data/MSDial/04-codeOutput/Thermo_results/volcano.csv"

read_file <- function(file_path, sep = ';') {
  # Check if file exists
  if (!file.exists(file_path)) {
    stop("File not found.")
  }
  
  # Extract file extension
  file_ext <- tools::file_ext(file_path)
  
  # Read file based on extension
  if (file_ext == "csv") {
    # Read CSV file
    data <- read.csv(file_path, stringsAsFactors = FALSE, sep=sep)
  } else if (file_ext == "xlsx") {
    # Read XLSX file
    data <- readxl::read_excel(file_path)
  } else {
    stop("Unsupported file format. Only CSV and XLSX files are supported.")
  }
  
  # Return the read data
  return(data)
}

filter_culumns <- function(data, column_names = NULL) {
  if (is.null(column_names)) {
    # Select columns starting with "M" if column_names is not specified
    filtered_columns <- data[, grep("^M", names(data)), drop = FALSE]
  } else {
    # Filter columns based on specified column names
    filtered_columns <- data[, column_names, drop = FALSE]
  }
  
  # Select columns starting with "M"
  m_columns <- grep("^M", names(data), value = TRUE)
  
  # Merge filtered columns with columns starting with "M"
  merged_data <- cbind(filtered_columns, data[, m_columns, drop = FALSE])
  
  return(merged_data)
}

save_as_csv <- function(dataframe, output_path, suffix) {
  # Extract directory path and filename without extension
  output_dir <- dirname(output_path)
  file_name <- paste0(tools::file_path_sans_ext(basename(output_path)), "_", suffix, ".csv")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Construct the full output file path
  output_file <- file.path(output_dir, file_name)
  
  # Save dataframe as CSV
  write.csv(dataframe, file = output_file, row.names = FALSE)
  
  cat("Dataframe saved as:", output_file, "\n")
}

## volcano plot function

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

generate_volcano_plot <- function(data, condition1_col, condition1_val, condition2_col, condition2_val) {
  # Filter metabolite data
  metabolite_data <- grep("^M", colnames(data), value = TRUE)
  
  # Extract condition 1 data
  condition1_data <- data %>%
    filter(!!sym(condition1_col) == condition1_val)
  
  # Extract condition 2 data
  condition2_data <- data %>%
    filter(!!sym(condition2_col) == condition2_val)
  
  # Check if conditions have enough observations
  if (nrow(condition1_data) == 0 || nrow(condition2_data) == 0) {
    stop("One or both conditions have no observations.")
  }
  
  # Calculate fold change
  fold_change <- sapply(metabolite_data, function(col) {
    mean_condition1 <- mean(condition1_data[[col]], na.rm = TRUE)
    mean_condition2 <- mean(condition2_data[[col]], na.rm = TRUE)
    fold_change <- mean_condition2 / mean_condition1
    return(fold_change)
  })
  
  fold_change_df <- data.frame(Metabolite = names(fold_change), Fold_Change = fold_change)
  
  # Calculate p-values
  p_values <- sapply(metabolite_data, function(col) {
    t_test_result <- t.test(condition1_data[[col]], condition2_data[[col]])
    return(t_test_result$p.value)
  })
  
  p_values_df <- data.frame(Metabolite = names(p_values), P_Value = p_values)
  
  # Merge fold change and p-value data frames
  volcano_df <- merge(p_values_df, fold_change_df, by = "Metabolite")
  volcano_df$log2FoldChange <- log2(volcano_df$Fold_Change)
  
  # Determine differentially expressed metabolites
  volcano_df$diffexpressed <- "NO"
  volcano_df$diffexpressed[volcano_df$log2FoldChange > 0.6 & volcano_df$P_Value < 0.05] <- "UP"
  volcano_df$diffexpressed[volcano_df$log2FoldChange < -0.6 & volcano_df$P_Value < 0.05] <- "DOWN"
  
  # Generate volcano plot
  volcano_plot <- ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(P_Value), col = diffexpressed, label = Metabolite)) +
    geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
    geom_point(size = 1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "orange") + # Add significance threshold line
    labs(color = 'Severe',x = "log2(FC)", y = "-log10(P-Value)") +
    theme_minimal()+
    theme(panel.background = element_rect(fill = "white", color = NA))+
    scale_color_manual(values = c("blue", "gray", "red"),
                       labels = c("Downregulated", "Not significant", "Upregulated"))
  
  return(list(volcano_plot = volcano_plot, volcano_data = volcano_df))
}





dataPOS <-read_file (InputcsvFilePOS)
dataNEG <-read_file (InputcsvFileNEG)


# PCA


# data numeric creation

X_dataPOS <-filter_culumns(dataPOS)
X_dataNEG <-filter_culumns(dataNEG)

# PCA processing
pca_resultPOS <- prcomp(X_dataPOS, scale = TRUE)
pca_resultNEG <- prcomp(X_dataNEG, scale = TRUE)

# PCA_POS dataframe creation
pca_dataPOS <- as.data.frame(pca_resultPOS$x)
pca_dataPOS$Etude <- dataPOS$Etude 
pca_dataPOS$sample_name <- dataPOS$sample_name
pca_dataPOS$Condition <- dataPOS$Condition
pca_dataPOS$ID <- dataPOS$ID
pca_dataPOS$Duration <- dataPOS$Duration
pca_dataPOS$SampleType <- dataPOS$SampleType
pca_dataPOS$batch <- dataPOS$batch
fpca_dataPOS = subset(pca_dataPOS, SampleType %in% c( 'sample'))
# Creer un dataframe pour les donnees de PCA_NEG
pca_dataNEG <- as.data.frame(pca_resultNEG$x)
pca_dataNEG$Etude <- dataNEG$Etude 
pca_dataNEG$sample_name <- dataNEG$sample_name
pca_dataNEG$Condition <- dataNEG$Condition
pca_dataNEG$ID <- dataNEG$ID
pca_dataNEG$Duration <- dataNEG$Duration
pca_dataNEG$SampleType <- dataNEG$SampleType
pca_dataNEG$batch <- dataNEG$batch
pca_dataNEG = subset(pca_dataNEG, SampleType %in% c('QC DIL', 'QC', 'blank', 'sample'))

# PCA processing
sample_colors <- c("sample" = "blue", "QC"="red", "QC DIL"="orange", "blank" = "black")
condition_colore <- c("T2D" = "blue", "Controle"="red", "Hirschsprung's disease"="orange")

# PCA pot creation
library(ggplot2)
library(ggrepel)
library(ggforce)
p_POS <- ggplot(fpca_dataPOS, aes(x = PC1, y = PC2, color = Condition, label = Duration)) +
  geom_point() +
  xlab("PC1") + ylab("PC2") +
  ggtitle("PCA POS") +
  scale_color_manual(values = condition_colore)  + 
  theme_minimal()+
  geom_text_repel(max.overlaps = Inf)+
  facet_wrap(~batch) #PCA decomposition 
print(p_POS)

p_NEG <- ggplot(pca_dataNEG, aes(x = PC1, y = PC2, color = Condition, label = ID)) +
  geom_point() +
  xlab("PC1") + ylab("PC2") +
  ggtitle("PCA NEG") +
  scale_color_manual(values = condition_colore) + 
  theme_minimal()+
  geom_text_repel(max.overlaps = Inf)
print(p_NEG)


# Volcano 


#POS
volcano_C_T_POS <- generate_volcano_plot(dataPOS, "Condition", "Controle", "Condition", "T2D")
volcano_C_H_POS <- generate_volcano_plot(dataPOS, "Condition", "Controle", "Condition", "Hirschsprung's disease")
volcano_T_H_POS <- generate_volcano_plot(dataPOS, "Condition", "T2D", "Condition", "Hirschsprung's disease")
volcano_C1_C2_POS <- generate_volcano_plot(dataPOS, "Duration", "C1", "Duration", "C2")
volcano_C1_T_POS <- generate_volcano_plot(dataPOS, "Duration", "C1", "Condition", "T2D")
volcano_C1_H_POS <- generate_volcano_plot(dataPOS, "Duration", "C1", "Condition", "Hirschsprung's disease")
volcano_C2_T_POS <- generate_volcano_plot(dataPOS, "Duration", "C2", "Condition", "T2D")
volcano_C2_H_POS <- generate_volcano_plot(dataPOS, "Duration", "C2", "Condition", "Hirschsprung's disease")

volcano_C_T_POS 
volcano_C_H_POS 
volcano_T_H_POS 
volcano_C1_C2_POS 
volcano_C1_T_POS 
volcano_C1_H_POS 
volcano_C2_T_POS 
volcano_C2_H_POS

save_as_csv(volcano_C_T_POS$volcano_data, OutputcsvFile, "C_T_POS")
save_as_csv(volcano_C_H_POS$volcano_data, OutputcsvFile, "C_H_POS")
save_as_csv(volcano_T_H_POS$volcano_data, OutputcsvFile, "T_H_POS")
save_as_csv(volcano_C1_C2_POS$volcano_data, OutputcsvFile, "C1_C2_POS")
save_as_csv(volcano_C1_T_POS$volcano_data, OutputcsvFile, "C1_T_POS")
save_as_csv(volcano_C1_H_POS$volcano_data, OutputcsvFile, "C1_H_POS")
save_as_csv(volcano_C2_T_POS$volcano_data, OutputcsvFile, "C2_T_POS")
save_as_csv(volcano_C2_H_POS$volcano_data, OutputcsvFile, "C2_H_POS")

# NEG
volcano_C_T_NEG <- generate_volcano_plot(dataNEG, "Condition", "Controle", "Condition", "T2D")
volcano_C_H_NEG <- generate_volcano_plot(dataNEG, "Condition", "Controle", "Condition", "Hirschsprung's disease")
volcano_T_H_NEG <- generate_volcano_plot(dataNEG, "Condition", "T2D", "Condition", "Hirschsprung's disease")
volcano_C1_C2_NEG <- generate_volcano_plot(dataNEG, "Duration", "C1", "Duration", "C2")
volcano_C1_T_NEG <- generate_volcano_plot(dataNEG, "Duration", "C1", "Condition", "T2D")
volcano_C1_H_NEG <- generate_volcano_plot(dataNEG, "Duration", "C1", "Condition", "Hirschsprung's disease")
volcano_C2_T_NEG <- generate_volcano_plot(dataNEG, "Duration", "C2", "Condition", "T2D")
volcano_C2_H_NEG <- generate_volcano_plot(dataNEG, "Duration", "C2", "Condition", "Hirschsprung's disease")

volcano_C_T_NEG 
volcano_C_H_NEG 
volcano_T_H_NEG 
volcano_C1_C2_NEG 
volcano_C1_T_NEG 
volcano_C1_H_NEG 
volcano_C2_T_NEG 
volcano_C2_H_NEG

save_as_csv(volcano_C_T_NEG$volcano_data, OutputcsvFile, "C_T_NEG")
save_as_csv(volcano_C_H_NEG$volcano_data, OutputcsvFile, "C_H_NEG")
save_as_csv(volcano_T_H_NEG$volcano_data, OutputcsvFile, "T_H_NEG")
save_as_csv(volcano_C1_C2_NEG$volcano_data, OutputcsvFile, "C1_C2_NEG")
save_as_csv(volcano_C1_T_NEG$volcano_data, OutputcsvFile, "C1_T_NEG")
save_as_csv(volcano_C1_H_NEG$volcano_data, OutputcsvFile, "C1_H_NEG")
save_as_csv(volcano_C2_T_NEG$volcano_data, OutputcsvFile, "C2_T_NEG")
save_as_csv(volcano_C2_H_NEG$volcano_data, OutputcsvFile, "C2_H_NEG")

