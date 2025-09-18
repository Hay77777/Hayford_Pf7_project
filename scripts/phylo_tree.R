# ===============================================================
# Script: phylogenetic_tree_from_snp_matrix.R
# Purpose: Generate a Neighbor-Joining phylogenetic tree 
#          from SNP genotype data.
# Author: [Hayford Osei Offei]
# Date: [18-09-2025]
# ===============================================================

# -----------------------------
# SETUP
# -----------------------------
setwd("E:/Malaria_Project/PLINK_analysis")  # Change to your working directory

# Load required libraries
library(dplyr)
library(tidyr)
library(ape)
library(tibble)

# -----------------------------
# STEP 1: LOAD AND PREPARE DATA
# -----------------------------
# Read the SNP data
snp_data <- read.delim(
  "PLINK_SNPs_summary_per_year_fixed.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Convert Genotypes to Binary (1 = SNP present, 0 = absent)
# Adjust logic here if your genotypes are coded differently (e.g., "0/0", "0/1", "1/1")
snp_data <- snp_data %>%
  mutate(Binary = ifelse(Genotype == "1/1", 1, 0))

# Create Unique Identifier combining Year and Sample_ID
# This ensures tree labels remain informative
snp_data <- snp_data %>%
  mutate(Unique_ID = paste(Year, Sample_ID, sep = "_"))  # e.g., "2010_ERR019549"

# -----------------------------
# STEP 2: CONVERT TO BINARY MATRIX
# -----------------------------
binary_snp_matrix <- snp_data %>%
  group_by(Unique_ID, SNP_Position) %>%          # Ensure unique combinations
  summarise(Binary = max(Binary), .groups = "drop") %>%  # Max = 1 if any SNP present
  pivot_wider(names_from = SNP_Position,
              values_from = Binary,
              values_fill = list(Binary = 0))

# Convert to matrix with Unique_ID as row names
binary_snp_matrix <- binary_snp_matrix %>%
  column_to_rownames(var = "Unique_ID") %>%
  as.matrix()

# -----------------------------
# STEP 3: COMPUTE DISTANCE MATRIX
# -----------------------------
# Use Jaccard distance (binary) for SNP presence/absence
dist_matrix <- dist(binary_snp_matrix, method = "binary")

# -----------------------------
# STEP 4: BUILD AND PLOT TREE
# -----------------------------
tree <- nj(dist_matrix)  # Neighbor-Joining method

# Retain full Unique_ID labels (Year_SampleID)
tree$tip.label <- rownames(binary_snp_matrix)

# Plot tree
plot(tree,
     main = "",
     cex = 0.35,          # Adjust label size
     label.offset = 0.02) # Move labels slightly away from branches

# -----------------------------
# STEP 5: SAVE OUTPUTS
# -----------------------------
# Save as PNG
png("Phylogenetic_Tree.png", width = 1200, height = 800, res = 150)
plot(tree, main = "", cex = 0.35, label.offset = 0.02)
dev.off()

# Save as PDF
pdf("Phylogenetic_Tree.pdf", width = 12, height = 8)
plot(tree, main = "", cex = 0.35, label.offset = 0.02)
dev.off()

cat("✅ Phylogenetic tree generated and saved as PNG + PDF.\n")
