# =====================================================================
#               AI and Biotechnology / Bioinformatics
# =====================================================================

# ---------------------------------------------------------------------
#              AI and Omics Research Internship (2025)


# Load packages
library(AnnotationDbi)
library(biomaRt)         # For probe-to-gene mapping (replacement for AnnotationDbi)
library(limma)           # Performs linear modeling and differential expression
library(dplyr)           # Simplifies data manipulation tasks
library(tibble)          
library(ggplot2)         # Used for plotting and visualization
library(pheatmap)        # Generates heatmaps for gene expression data

# -------------------------------------------------------------
#### Load Preprocessed Data ####
# -------------------------------------------------------------

load("E-MTAB-9990.Rdata")

# Check annotation slot of your dataset
annotation(raw_data)

raw_data

# -------------------------------------------------------------
#### Probe-to-Gene Mapping using biomaRt ####
# -------------------------------------------------------------

# For PrimeView platform, we use biomaRt instead of AnnotationDbi
# because there is no primeview.db annotation package available

cat("=== CONNECTING TO ENSEMBL DATABASE ===\n")
# Connect to Ensembl biomaRt
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "hsapiens_gene_ensembl",
                      mirror = "useast")

# Check available filters
filters <- listFilters(ensembl)
primeview_filters <- grep("primeview", filters$name, value = TRUE, ignore.case = TRUE)
cat("Available PrimeView filters:\n")
print(primeview_filters)
cat("\n")

# -------------------------------------------------------------
# Extract probe IDs from processed microarray data
# -------------------------------------------------------------
probe_ids <- rownames(processed_data)
cat("Total number of probes:", length(probe_ids), "\n\n")

# Map probe IDs to gene symbols using biomaRt
cat("Mapping probes to genes (this may take 1-2 minutes)...\n")

# Process in batches to avoid timeout
batch_size <- 5000
n_batches <- ceiling(length(probe_ids) / batch_size)
gene_mapping_list <- list()

for (i in 1:n_batches) {
  start_idx <- (i - 1) * batch_size + 1
  end_idx <- min(i * batch_size, length(probe_ids))
  batch_probes <- probe_ids[start_idx:end_idx]
  
  cat("Processing batch", i, "of", n_batches, "...\n")
  
  tryCatch({
    batch_mapping <- getBM(
      attributes = c('affy_primeview', 'hgnc_symbol'),
      filters = 'affy_primeview',
      values = batch_probes,
      mart = ensembl
    )
    gene_mapping_list[[i]] <- batch_mapping
  }, error = function(e) {
    cat("Error in batch", i, ":", e$message, "\n")
  })
}

# Combine all batches
gene_mapping <- do.call(rbind, gene_mapping_list)
colnames(gene_mapping) <- c("PROBEID", "SYMBOL")

cat("\nMapping complete!\n")
cat("Probes successfully mapped:", nrow(gene_mapping), "\n")

# Create a complete mapping with NAs for unmapped probes
gene_symbols_df <- data.frame(
  PROBEID = probe_ids,
  stringsAsFactors = FALSE
) %>%
  left_join(gene_mapping, by = "PROBEID") %>%
  mutate(SYMBOL = ifelse(SYMBOL == "" | is.na(SYMBOL), NA, SYMBOL))

# Convert to named vector (same format as mapIds output)
gene_symbols <- setNames(gene_symbols_df$SYMBOL, gene_symbols_df$PROBEID)

# Convert mapping to a data frame and rename columns
gene_map_df <- gene_symbols_df

cat("Probes with gene symbols:", sum(!is.na(gene_map_df$SYMBOL)), "\n")
cat("Probes without gene symbols:", sum(is.na(gene_map_df$SYMBOL)), "\n\n")

# -------------------------------------------------------------
# Handle multiple probes mapping to a single gene
# -------------------------------------------------------------
# Several strategies exist:
# 1. Retain probe with highest expression or variance
# 2. Average or summarize probe signals
# 3. Remove duplicate probes to maintain one row per gene

# Summarize number of probes per gene symbol
duplicate_summary <- gene_map_df %>%
  filter(!is.na(SYMBOL)) %>%
  group_by(SYMBOL) %>%
  summarise(probes_per_gene = n()) %>%
  arrange(desc(probes_per_gene))

# Identify genes associated with multiple probes
duplicate_genes <- duplicate_summary %>%
  filter(probes_per_gene > 1)

cat("=== DUPLICATE PROBE ANALYSIS ===\n")
cat("Total unique genes:", nrow(duplicate_summary), "\n")
cat("Genes with multiple probes:", nrow(duplicate_genes), "\n")
cat("Total duplicate probes:", sum(duplicate_genes$probes_per_gene), "\n")
cat("Max probes per gene:", max(duplicate_summary$probes_per_gene), "\n\n")

cat("Top 10 genes with most probes:\n")
print(head(duplicate_genes, 10))
cat("\n")

# -------------------------------------------------------------
# Merge annotation mapping with expression data
# -------------------------------------------------------------

# Verify if probe IDs in mapping correspond to expression data
all(gene_map_df$PROBEID == row.names(processed_data))

# Merge annotation (SYMBOL) with expression matrix
processed_data_df <- processed_data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::mutate(SYMBOL = gene_symbols[PROBEID]) %>%
  dplyr::relocate(SYMBOL, .after = PROBEID)

# Remove probes without valid gene symbol annotation
processed_data_df <- processed_data_df %>%
  dplyr::filter(!is.na(SYMBOL))

cat("Probes remaining after filtering NAs:", nrow(processed_data_df), "\n")

# Select only numeric expression columns
expr_only <- processed_data_df %>%
  dplyr::select(-PROBEID, -SYMBOL)

# -------------------------------------------------------------
# Collapse multiple probes per gene using average expression
# -------------------------------------------------------------
# limma::avereps() computes the average for probes representing the same gene
averaged_data <- limma::avereps(expr_only, ID = processed_data_df$SYMBOL)

cat("Unique genes after averaging:", nrow(averaged_data), "\n\n")

# Convert averaged expression data to matrix format
data <- as.data.frame(averaged_data)
data <- data.matrix(data)
str(data)        # Structure check
is.numeric(data) # Confirm numeric matrix

# -------------------------------------------------------------
#### Differential Gene Expression Analysis ####
# -------------------------------------------------------------

# Define sample groups based on phenotype data
# For E-MTAB-9990: tumor vs normal tissue adjacent to tumor

# Check phenotype column names
cat("=== PHENOTYPE DATA ===\n")
cat("Column names:\n")
print(colnames(phenotype_data))
cat("\n")

# Identify the correct column for sample groups
group_col <- "Characteristics.sampling.site."

cat("Using column:", group_col, "\n")
cat("Unique values:\n")
print(unique(phenotype_data[[group_col]]))
cat("\n")

groups <- factor(phenotype_data[[group_col]],
                 levels = c("normal tissue adjacent to tumor", "tumor"),
                 labels = c("normal", "cancer"))

class(groups)
levels(groups)

cat("Sample distribution:\n")
print(table(groups))
cat("\n")

# -------------------------------------------------------------
# Create design matrix for linear modeling
# -------------------------------------------------------------
# Using no intercept (~0 + groups) allows each group to have its own coefficient
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

cat("Design matrix dimensions:", dim(design), "\n")
cat("Data matrix dimensions:", dim(data), "\n\n")

# Fit linear model to expression data
fit_1 <- lmFit(data, design)

# Define contrast to compare cancer vs normal samples
contrast_matrix <- makeContrasts(cancer_vs_normal = cancer - normal,
                                 levels = design)

# Apply contrasts and compute moderated statistics
fit_contrast <- contrasts.fit(fit_1, contrast_matrix)

fit_2 <- eBayes(fit_contrast)

# -------------------------------------------------------------
# Extract list of differentially expressed genes (DEGs)
# -------------------------------------------------------------
deg_results <- topTable(fit_2,
                        coef = "cancer_vs_normal",  # Specify contrast of interest
                        number = Inf,               # Return all genes
                        adjust.method = "BH")       # Benjamini-Hochberg correction

# -------------------------------------------------------------
# Classify DEGs into Upregulated, Downregulated, or Not Significant
# -------------------------------------------------------------
deg_results$threshold <- as.factor(ifelse(
  deg_results$adj.P.Val < 0.05 & deg_results$logFC > 1, "Upregulated",
  ifelse(deg_results$adj.P.Val < 0.05 & deg_results$logFC < -1, "Downregulated",
         "Not Significant")
))

# Subset genes by regulation direction
upregulated <- subset(deg_results, threshold == "Upregulated")
downregulated <- subset(deg_results, threshold == "Downregulated")

# Combine both sets of DEGs
deg_updown <- rbind(upregulated, downregulated)

cat("=== DIFFERENTIAL EXPRESSION RESULTS ===\n")
cat("Total genes analyzed:", nrow(deg_results), "\n")
cat("Upregulated genes (logFC > 1, adj.P.Val < 0.05):", nrow(upregulated), "\n")
cat("Downregulated genes (logFC < -1, adj.P.Val < 0.05):", nrow(downregulated), "\n")
cat("Total significant DEGs:", nrow(deg_updown), "\n\n")

# Create Results directory
dir.create("Results", showWarnings = FALSE)
dir.create("Result_Plots", showWarnings = FALSE)

write.csv(deg_results, file = "Results/DEGs_Results.csv", row.names = TRUE)
write.csv(upregulated, file = "Results/Upregulated_DEGs.csv", row.names = TRUE)
write.csv(downregulated, file = "Results/Downregulated_DEGs.csv", row.names = TRUE)
write.csv(deg_updown, file = "Results/Updown_DEGs.csv", row.names = TRUE)

cat("Results saved to CSV files\n\n")

# -------------------------------------------------------------
#### Data Visualization ####
# -------------------------------------------------------------

# -------------------------------------------------------------
# Volcano Plot: visualizes DEGs by logFC and adjusted p-values
# -------------------------------------------------------------
# Note: x-axis = log2 fold change, y-axis = -log10 adjusted p-value

# Save volcano plot as PNG
png("Result_Plots/volcano_plot1.png", width = 2000, height = 1500, res = 300)

ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       subtitle = "Gastric Cancer Tumor vs Normal Adjacent Tissue",
       x = "log2 Fold Change",
       y = "-log10(Adjusted P-value)",
       color = "Regulation") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black")

dev.off()

cat("Volcano plot saved\n")

# -------------------------------------------------------------
# Heatmap of Top 25 Differentially Expressed Genes
# -------------------------------------------------------------

# Select top 25 genes with smallest adjusted p-values
n_top <- min(25, nrow(deg_updown))
top_genes <- head(rownames(deg_updown[order(deg_updown$adj.P.Val), ]), n_top)

# Subset averaged expression matrix for selected genes
heatmap_data <- data[top_genes, ]

# Generate unique column names per sample group for display
group_char <- as.character(groups)
heatmap_names <- ave(group_char, group_char, FUN = function(x) paste0(x, "_", seq_along(x)))

# Assign formatted names to heatmap columns
colnames(heatmap_data) <- heatmap_names

# Create annotation for heatmap
annotation_col <- data.frame(
  Group = groups,
  row.names = colnames(heatmap_data)
)

# Save heatmap as PNG
png("Result_Plots/heatmap_top25_DEGs1.png", width = 2000, height = 1500, res = 300)

# Generate heatmap with row scaling
pheatmap(
  heatmap_data,
  scale = "row",                     # Scale by row for better visualization
  cluster_rows = TRUE,               # Cluster genes
  cluster_cols = TRUE,               # Cluster samples
  show_rownames = TRUE,              # Display gene names
  show_colnames = TRUE,              # Display sample labels
  annotation_col = annotation_col,   # Add group annotation
  color = colorRampPalette(c("blue", "white", "red"))(100),
  fontsize_row = 8,
  fontsize_col = 8,
  main = paste0("Top ", n_top, " Differentially Expressed Genes")
)

dev.off()

cat("Heatmap saved\n\n")

save.image("class3c assignment.Rdata")

# =====================================================================
# END OF SCRIPT
# =====================================================================