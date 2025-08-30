# --------------------------
# Assignment 2
# --------------------------
# In this assignment you will work with the results of differential gene expression (DGE) analysis. 
#The analysis produces two key measures for each gene:

# log2FoldChange (log2FC): 
# Indicates the magnitude and direction of change in gene expression. 
# Positive values suggest higher expression(upregulated gene) in the experimental condition compared to control. 
# Negative values suggest lower expression (downregulated gene). 
# The absolute value reflects the strength of the change.

# Adjusted p-value (padj): 
# Represents the statistical significance of the observed difference, corrected for multiple testing. 
# A smaller value indicates stronger evidence that the observed difference is not due to chance.

# Write a function classify_gene() 

# that takes:
#   - logFC (log2FoldChange)
#   - padj  (adjusted p-value)

# and returns:
#   - "Upregulated" if log2FC > 1 and padj < 0.05
#   - "Downregulated" if log2FC < -1 and padj < 0.05
#   - "Not_Significant" otherwise

# function classify_gene

Classify_gene <- function(logFC, padj) {
  ifelse(logFC > 1 & padj < 0.05, "Upregulated",
         ifelse(logFC < -1 & padj < 0.05, "Downregulated", "Not_Significant"))
}


# Then:
#   - Apply it in a for-loop to process both datasets (DEGs_data_1.csv, DEGs_data_2.csv)
#   - Replace missing padj values with 1
#   - Add a new column 'status'
#   - Save processed files into Results folder
#   - Print summary counts of significant, upregulated, and downregulated genes
#   - Use table() for summaries

# Defining the input folder and the output folder 
input_dir_1 <- "raw_data"
output_dir_1 <- "results"

# listing the files to process
files_to_process_1 <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv")

# Preparing an empty list in R to store results for later use.
result_list_1 <- list()

# for loop to process the datasets
for (file_names_1 in files_to_process_1) {
  cat("\nProcessing:", file_names_1, "\n")
  
  input_file_path_1 <- file.path(input_dir, file_names_1)
  
  #import the datasets
  data_1 <- read.csv(input_file_path_1, header = TRUE)
  cat("File imported. Checking for missing values... \n")
  
  # handling missing values in padj column
  if("padj" %in% names(data_1)){
    missing_count <- sum(is.na(data_1$padj))
    
    cat("Missing values in 'padj':", missing_count, "\n")
    data_1$padj[is.na(data_1$padj)] <- mean(data_1$padj, na.rm = TRUE)
  }
  
  # Checking gene classification
  data_1$status <-Classify_gene(data_1$logFC, data_1$padj)
  cat("Gene status is sucessfully identified. \n")
  
  # saving results in R
  result_list_1[[file_names_1]] <- data_1
  
  # saving results in Results folder
  output_file_path_1 <-file.path(output_dir_1, paste0("Status results", file_names_1))
  write.csv(data_1, output_file_path_1, row.names = FALSE)
  cat("Results saved to:", output_file_path_1, "\n")
  
}

#  Print summary counts of significant, upregulated, and downregulated genes
#   Use table() for summaries
print(table(data_1$status))

