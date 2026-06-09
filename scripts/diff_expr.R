#!/usr/bin/env Rscript
# diff_expr.R
# Conducts differential expression analysis using DESeq2 on Salmon outputs.

cat("=== Starting Differential Expression Analysis in R ===\n")

# Install pacman if not available to easily handle library loading
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

# Load required libraries
library(DESeq2)
library(ggplot2)

# Establish directories
scripts_dir <- getwd() # script runs with working dir in rnaseq_pipeline/
quants_dir <- file.path(scripts_dir, "quants")
output_dir <- scripts_dir

# 1. Load Salmon quantification files
samples <- c("Control1", "Control2", "Control3", "Treated1", "Treated2", "Treated3")
quant_files <- file.path(quants_dir, paste0(samples, "_quant"), "quant.sf")

# Verify all files exist
for (f in quant_files) {
  if (!file.exists(f)) {
    stop(paste("Error: Quantification file not found at", f))
  }
}

cat("Loading Salmon quantification files...\n")

# Read the first file to get gene list
first_quant <- read.table(quant_files[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
genes <- first_quant$Name
num_genes <- length(genes)

# Initialize counts matrix
count_matrix <- matrix(0, nrow = num_genes, ncol = length(samples))
rownames(count_matrix) <- genes
colnames(count_matrix) <- samples

# Fill counts matrix with estimated NumReads from Salmon
for (i in 1:length(samples)) {
  quant_data <- read.table(quant_files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # Match names to ensure alignment
  count_matrix[, i] <- quant_data$NumReads[match(genes, quant_data$Name)]
}

# DESeq2 requires integer counts (Salmon outputs decimals for estimated reads)
count_matrix <- round(count_matrix)

cat("Counts matrix loaded. Dimensions:", paste(dim(count_matrix), collapse = "x"), "\n")

# 2. Setup experimental design Metadata (colData)
col_data <- data.frame(
  row.names = samples,
  condition = factor(c("control", "control", "control", "treated", "treated", "treated"),
                     levels = c("control", "treated"))
)

cat("Experimental Design:\n")
print(col_data)

# 3. Create DESeq2 Dataset and run analysis
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = col_data,
                              design = ~ condition)

cat("Running DESeq2...\n")
dds <- DESeq(dds)
res <- results(dds)

# Sort results by adjusted p-value
res_ordered <- res[order(res$padj), ]
res_df <- as.data.frame(res_ordered)

# Save results to CSV
output_csv <- file.path(output_dir, "diff_expr_results.csv")
write.csv(res_df, file = output_csv)
cat("Saved differential expression results to:", output_csv, "\n")

# 4. Generate Volcano Plot
cat("Generating Volcano Plot...\n")
res_df$gene <- rownames(res_df)

# Prevent Inf values on Y-axis by capping p-value at 1e-250
res_df$pvalue[res_df$pvalue == 0] <- 1e-250

# Define significance thresholds
padj_cutoff <- 0.05
lfc_cutoff <- 1.0

# Classify genes for coloring
res_df$Significance <- "Not Significant"
res_df$Significance[res_df$padj < padj_cutoff & res_df$log2FoldChange > lfc_cutoff] <- "Upregulated"
res_df$Significance[res_df$padj < padj_cutoff & res_df$log2FoldChange < -lfc_cutoff] <- "Downregulated"
res_df$Significance <- factor(res_df$Significance, levels = c("Not Significant", "Upregulated", "Downregulated"))

# Plot
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = Significance)) +
  geom_point(alpha = 0.8, size = 2.5) +
  scale_color_manual(values = c("Not Significant" = "grey60", "Upregulated" = "firebrick3", "Downregulated" = "dodgerblue3")) +
  theme_bw(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey40") +
  labs(title = "Volcano Plot: Treated vs. Control",
       x = "Log2 Fold Change",
       y = "-Log10 P-value")

# Label the top differentially expressed genes
top_genes <- head(res_df[res_df$Significance != "Not Significant", ], 6)
if (nrow(top_genes) > 0) {
  volcano_plot <- volcano_plot +
    geom_text(data = top_genes, aes(label = gene), vjust = -0.6, hjust = 0.5, color = "black", size = 3.5, fontface = "bold", show.legend = FALSE)
}

# Save plot to PNG
output_png <- file.path(output_dir, "volcano_plot.png")
ggsave(output_png, plot = volcano_plot, width = 8, height = 6, dpi = 300)
cat("Saved Volcano Plot to:", output_png, "\n")

cat("=== Differential Expression Completed Successfully ===\n")
