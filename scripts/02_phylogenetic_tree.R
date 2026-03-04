# ===============================================================
# Script:  02_phylogenetic_tree.R
# Purpose: Generate a Neighbor-Joining phylogenetic tree from
#          SNP genotype data using binary presence/absence matrix
# Study:   Longitudinal Genomic Surveillance of Antimalarial
#          Drug Resistance-Associated SNPs, Upper East Ghana
# Author:  Hayford Osei Offei
# Date:    2025
# ===============================================================
#
# INPUT FILES REQUIRED:
#   - PLINK_SNPs_summary_per_year_fixed.tsv
#       Columns: Year, Sample_ID, CHROM, SNP_Position, REF, ALT,
#                Gene_ID, Amino_Acid_Change, Genotype, Count
#
# OUTPUT FILES:
#   - results/Phylogenetic_Tree.png
#   - results/Phylogenetic_Tree.pdf
#
# DEPENDENCIES:
#   - R packages: dplyr, tidyr, tibble, ape, ggplot2
# ===============================================================


# ---------------------------------------------------------------
# SETUP
# ---------------------------------------------------------------
setwd("E:/Malaria_Project/PLINK_analysis")  # Update path as needed

dir.create("results", showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ape)
  library(ggplot2)
})


# ---------------------------------------------------------------
# STEP 1: LOAD AND PREPARE SNP DATA
# ---------------------------------------------------------------
snp_data <- read.delim(
  "PLINK_SNPs_summary_per_year_fixed.tsv",
  header          = TRUE,
  sep             = "\t",
  stringsAsFactors = FALSE
)

cat("Loaded SNP data:", nrow(snp_data), "rows,",
    length(unique(snp_data$Sample_ID)), "samples\n")

# Convert genotype to binary (1 = SNP present, 0 = absent/reference)
# Genotype "1/1" = homozygous alternate (SNP present)
# Genotype "0/0" or other = reference/missing (SNP absent)
snp_data <- snp_data %>%
  mutate(
    Binary    = ifelse(Genotype == "1/1", 1, 0),
    Unique_ID = paste(Year, Sample_ID, sep = "_")  # e.g. "2010_ERR019549"
  )


# ---------------------------------------------------------------
# STEP 2: BUILD BINARY SNP MATRIX (samples x SNP positions)
# ---------------------------------------------------------------
# Each row = one sample, each column = one SNP genomic position
# Cell value = 1 if SNP present in that sample, 0 if absent
binary_matrix <- snp_data %>%
  group_by(Unique_ID, SNP_Position) %>%
  summarise(Binary = max(Binary), .groups = "drop") %>%
  pivot_wider(
    names_from   = SNP_Position,
    values_from  = Binary,
    values_fill  = list(Binary = 0)
  ) %>%
  column_to_rownames(var = "Unique_ID") %>%
  as.matrix()

cat("Binary SNP matrix dimensions:", nrow(binary_matrix),
    "samples x", ncol(binary_matrix), "SNP positions\n")


# ---------------------------------------------------------------
# STEP 3: COMPUTE GENETIC DISTANCE MATRIX
# ---------------------------------------------------------------
# Jaccard (binary) distance: appropriate for binary
# presence/absence SNP data
dist_matrix <- dist(binary_matrix, method = "binary")


# ---------------------------------------------------------------
# STEP 4: BUILD NEIGHBOR-JOINING TREE
# ---------------------------------------------------------------
tree <- nj(dist_matrix)

# Ensure tip labels match the binary matrix row names (Year_SampleID)
tree$tip.label <- rownames(binary_matrix)

# Extract year from tip label for colour annotation
tip_years <- as.integer(sub("_.*", "", tree$tip.label))
year_palette <- colorRampPalette(c("#2166AC", "#74ADD1", "#ABD9E9",
                                    "#FEE090", "#FDAE61", "#F46D43",
                                    "#D73027", "#A50026", "#67001F"))(
  length(unique(tip_years))
)
tip_colors <- year_palette[as.numeric(factor(tip_years))]


# ---------------------------------------------------------------
# STEP 5: SAVE OUTPUTS
# ---------------------------------------------------------------

# --- PNG output ---
png("results/Phylogenetic_Tree.png",
    width = 1400, height = 1000, res = 150)

plot(tree,
     type         = "phylogram",
     cex          = 0.35,
     label.offset = 0.005,
     tip.color    = tip_colors,
     edge.color   = "gray30",
     edge.width   = 0.8,
     main         = "")

# Add year legend
legend("bottomleft",
       legend = sort(unique(tip_years)),
       fill   = year_palette,
       title  = "Collection Year",
       cex    = 0.6,
       bty    = "n")

dev.off()
cat("Phylogenetic tree (PNG) saved to results/Phylogenetic_Tree.png\n")


# --- PDF output (vector, publication quality) ---
pdf("results/Phylogenetic_Tree.pdf", width = 14, height = 10)

plot(tree,
     type         = "phylogram",
     cex          = 0.35,
     label.offset = 0.005,
     tip.color    = tip_colors,
     edge.color   = "gray30",
     edge.width   = 0.8,
     main         = "")

legend("bottomleft",
       legend = sort(unique(tip_years)),
       fill   = year_palette,
       title  = "Collection Year",
       cex    = 0.6,
       bty    = "n")

dev.off()
cat("Phylogenetic tree (PDF) saved to results/Phylogenetic_Tree.pdf\n")


# ---------------------------------------------------------------
# STEP 6: BASIC TREE STATISTICS
# ---------------------------------------------------------------
cat("\nTree summary:\n")
cat("  Number of tips  :", length(tree$tip.label), "\n")
cat("  Number of nodes :", tree$Nnode, "\n")
cat("  Is rooted       :", is.rooted(tree), "\n")

# Per-year tip count
tip_df <- data.frame(
  Tip  = tree$tip.label,
  Year = tip_years
)
cat("\nSamples per year in tree:\n")
print(table(tip_df$Year))
