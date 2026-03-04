# ===============================================================
# Script:  01_plink_pca_analysis.R
# Purpose: PLINK-based SNP QC, PCA, hierarchical clustering,
#          and statistical analysis of P. falciparum WGS data
# Study:   Longitudinal Genomic Surveillance of Antimalarial
#          Drug Resistance-Associated SNPs, Upper East Ghana
# Author:  Hayford Osei Offei
# Date:    2025
# ===============================================================
#
# INPUT FILES REQUIRED:
#   - Merged VCF file of all 77 QC-passed samples (merged.vcf.gz)
#   - PLINK binary files generated in Steps 1-3 below
#
# OUTPUT FILES:
#   - pf_snp_pca.eigenvec       (PCA coordinates)
#   - pf_snp_pca.eigenval       (variance explained per PC)
#   - results/PCA_plot.png
#   - results/Hierarchical_Clustering_Dendrogram.png
#
# DEPENDENCIES:
#   - PLINK v1.9 (command line, run before this R script)
#   - R packages: ggplot2, ggrepel, ggdendro, cluster, ape
# ===============================================================


# ---------------------------------------------------------------
# PRE-STEP: PLINK COMMANDS (run in terminal before this R script)
# ---------------------------------------------------------------
# These bash commands were run in WSL2 Ubuntu on an Intel Core i7
# system with 16 GB RAM. Run them in your terminal, NOT in R.
#
# Step A: Convert PED/MAP files to PLINK binary format
#         Input: pf_snp.ped + pf_snp.map
#         (69 variants, 71 samples)
#         Note: Variant 37 was triallelic; rarest allele set missing
#
#   plink --file pf_snp \
#         --make-bed \
#         --out pf_snp_binary
#
# Step B: Impute missing genotypes using the reference allele (A2)
#         Total genotyping rate was 0.1357 before imputation
#
#   plink --bfile pf_snp_binary \
#         --fill-missing-a2 \
#         --make-bed \
#         --out pf_snp_imputed
#
# Step C: Run PCA on imputed data (10 principal components)
#         Final dataset: 69 variants, 71 samples
#         Output: pf_snp_pca.eigenval + pf_snp_pca.eigenvec
#
#   plink --bfile pf_snp_imputed \
#         --pca 10 \
#         --out pf_snp_pca
#
# After running the above, you will have:
#   pf_snp_pca.eigenvec   <- used in STEP 1 below
#   pf_snp_pca.eigenval   <- used in STEP 3 below
#
# NOTE ON SAMPLE NUMBERS:
#   90 samples downloaded from ENA
#   82 samples passed FastQC/Trimmomatic QC
#   71 samples retained after PLINK processing
#   Per-year breakdown of retained samples:
#     2010: 8 | 2011: 6 | 2012: 6 | 2013: 9 | 2014: 8
#     2015: 10 | 2016: 10 | 2017: 4 | 2018: 10
# ---------------------------------------------------------------


# ---------------------------------------------------------------
# SETUP
# ---------------------------------------------------------------
setwd("E:/Malaria_Project/PLINK_analysis")  # Update path as needed

# Create results directory if it doesn't exist
dir.create("results", showWarnings = FALSE)

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(ggdendro)
  library(cluster)
  library(ape)
})


# ---------------------------------------------------------------
# STEP 1: LOAD PCA DATA
# ---------------------------------------------------------------
# Load eigenvec file produced by PLINK --pca
pca_data <- read.table("pf_snp_pca.eigenvec", header = FALSE)
colnames(pca_data) <- c("FID", "Sample_ID", paste0("PC", 1:(ncol(pca_data) - 2)))

cat("PCA data loaded:", nrow(pca_data), "samples,",
    ncol(pca_data) - 2, "PCs\n")
print(head(pca_data[, 1:6]))


# ---------------------------------------------------------------
# STEP 2: LOAD SAMPLE METADATA (year labels)
# ---------------------------------------------------------------
# Add sampling year to each sample for colouring the PCA plot.
# This maps ERR accession IDs to their collection year.
year_map <- data.frame(
  Sample_ID = c(
    # 2010
    "ERR019549","ERR039235","ERR036147","ERR035362","ERR045604",
    "ERR205933","ERR211562","ERR045620","ERR045607","ERR045609",
    # 2011
    "ERR039979","ERR042698","ERR042701","ERR042719","ERR042722",
    "ERR042727","ERR063567","ERR114406","ERR114409","ERR211462",
    # 2012
    "ERR376187","ERR376188","ERR376201","ERR376202","ERR376203",
    "ERR376212","ERR376213","ERR376214","ERR376216","ERR403196",
    # 2013
    "ERR450116","ERR450108","ERR450099","ERR450045","ERR586212",
    "ERR586179","ERR450104","ERR450112","ERR450044","ERR450041",
    # 2014
    "ERR3504126","ERR3523699","ERR3546199","ERR3546151","ERR3546188",
    "ERR3594524","ERR3594529","ERR3523718","ERR3523721","ERR3594576",
    # 2015
    "ERR3486201","ERR3523571","ERR3523563","ERR3546220","ERR3523635",
    "ERR3523614","ERR1214163","ERR1214188","ERR1214209","ERR1214148",
    # 2016
    "ERR2532525","ERR2542018","ERR2541936","ERR2532503","ERR2532517",
    "ERR2532524","ERR2532570","ERR2541897","ERR2541942","ERR2532534",
    # 2017
    "ERR3486300","ERR2891366","ERR2890914","ERR3608767","ERR3546081",
    "ERR2375226","ERR2506008","ERR3504084","ERR3523658","ERR3504120",
    # 2018
    "ERR3608800","ERR3608865","ERR3584487","ERR3594738","ERR3608785",
    "ERR3594774","ERR3594757","ERR3594689","ERR3594679","ERR3594661"
  ),
  Year = c(
    rep(2010, 10), rep(2011, 10), rep(2012, 10), rep(2013, 10),
    rep(2014, 10), rep(2015, 10), rep(2016, 10), rep(2017, 10),
    rep(2018, 10)
  )
)

# Merge year into PCA data
pca_data <- merge(pca_data, year_map, by = "Sample_ID", all.x = TRUE)
pca_data$Year <- as.factor(pca_data$Year)


# ---------------------------------------------------------------
# STEP 3: VARIANCE EXPLAINED BY EACH PC
# ---------------------------------------------------------------
eigenvalues      <- scan("pf_snp_pca.eigenval")
variance_explained <- round((eigenvalues / sum(eigenvalues)) * 100, 2)

cat("PC1 explains:", variance_explained[1], "%\n")
cat("PC2 explains:", variance_explained[2], "%\n")

pc1_label <- paste0("PC1 (", variance_explained[1], "%)")
pc2_label <- paste0("PC2 (", variance_explained[2], "%)")


# ---------------------------------------------------------------
# STEP 4: IDENTIFY OUTLIERS (top 5% by Euclidean distance)
# ---------------------------------------------------------------
pca_data$distance <- sqrt(pca_data$PC1^2 + pca_data$PC2^2)
threshold         <- quantile(pca_data$distance, 0.95, na.rm = TRUE)
pca_data$outlier  <- pca_data$distance > threshold

cat("Outlier threshold (95th percentile):", round(threshold, 4), "\n")
cat("Number of outliers:", sum(pca_data$outlier), "\n")


# ---------------------------------------------------------------
# STEP 5: PCA SCATTER PLOT (coloured by year)
# ---------------------------------------------------------------
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2,
                                  color = Year,
                                  shape = outlier)) +
  geom_jitter(width = 0.02, height = 0.02, size = 2, alpha = 0.75) +
  scale_shape_manual(values = c(16, 17),
                     labels = c("Main cluster", "Outlier")) +
  geom_text_repel(
    data    = subset(pca_data, outlier),
    aes(label = Sample_ID),
    size    = 2.5,
    color   = "black",
    segment.color = "grey60"
  ) +
  labs(
    x      = pc1_label,
    y      = pc2_label,
    color  = "Collection Year",
    shape  = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position   = "right",
    panel.grid.major  = element_line(color = "gray90"),
    panel.grid.minor  = element_blank()
  )

ggsave("results/PCA_plot.png", plot = pca_plot,
       width = 7, height = 5, dpi = 300)
cat("PCA plot saved to results/PCA_plot.png\n")


# ---------------------------------------------------------------
# STEP 6: HIERARCHICAL CLUSTERING ON PCA COORDINATES
# ---------------------------------------------------------------
# Use first 10 PCs for clustering
pc_cols       <- grep("^PC", colnames(pca_data), value = TRUE)[1:10]
dist_matrix   <- dist(pca_data[, pc_cols], method = "euclidean")
hc            <- hclust(dist_matrix, method = "ward.D2")

# Base R dendrogram (high-resolution)
png("results/Hierarchical_Clustering_Dendrogram.png",
    width = 10, height = 6, units = "in", res = 300)
plot(hc,
     labels = pca_data$Sample_ID,
     main   = "",
     xlab   = "",
     sub    = "",
     cex    = 0.35)
dev.off()
cat("Dendrogram saved to results/Hierarchical_Clustering_Dendrogram.png\n")

# ggplot2 dendrogram
dendro_data    <- ggdendro::dendro_data(as.dendrogram(hc))
dendrogram_plot <- ggplot() +
  geom_segment(data = dendro_data$segments,
               aes(x = x, y = y, xend = xend, yend = yend),
               color = "steelblue", linewidth = 0.5) +
  theme_minimal() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title   = element_blank()
  )

ggsave("results/Dendrogram_ggplot.png", plot = dendrogram_plot,
       width = 10, height = 4, dpi = 300)


# ---------------------------------------------------------------
# STEP 7: STATISTICAL ANALYSIS
# ---------------------------------------------------------------
# ANOVA: test whether outlier samples differ significantly in
# Euclidean distance from the main cluster centroid
pca_data$group <- ifelse(pca_data$outlier, "Outlier", "Main Cluster")
anova_result   <- aov(distance ~ group, data = pca_data)
cat("\nANOVA Result (outlier vs main cluster):\n")
print(summary(anova_result))

# Cophenetic correlation coefficient (measures how well the
# hierarchical clustering preserves pairwise distances)
ccc <- cor(dist_matrix, cophenetic(hc))
cat("\nCophenetic Correlation Coefficient:", round(ccc, 4), "\n")
