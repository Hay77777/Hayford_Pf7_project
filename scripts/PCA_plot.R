# ===============================================================
# Script: pca_and_clustering_analysis.R
# Purpose: Perform PCA visualization, hierarchical clustering, 
#          dendrogram plotting, and statistical analysis of SNP data.
# Author: [Hayford Osei Offei]
# Date: [18-09-2025]
# ===============================================================

# -----------------------------
# SETUP
# -----------------------------
setwd("E:/Malaria_Project/PLINK_analysis")  # Change to working directory

# Load required libraries
suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
  library(ggrepel)
  library(ggdendro)
  library(cluster)
})

# -----------------------------
# STEP 1: LOAD PCA DATA
# -----------------------------
pca_data <- read.table("pf_snp_pca.eigenvec", header = FALSE)
colnames(pca_data) <- c("Sample", "Individual", paste0("PC", 1:(ncol(pca_data) - 2)))

# Preview
cat("✅ PCA data loaded:\n")
print(head(pca_data))

# -----------------------------
# STEP 2: PCA PLOT (PC1 vs PC2)
# -----------------------------
# Compute Euclidean distance from (0,0) for PC1 & PC2
pca_data$distance <- sqrt(pca_data$PC1^2 + pca_data$PC2^2)

# Identify outliers (top 5% farthest)
threshold <- quantile(pca_data$distance, 0.95, na.rm = TRUE)
pca_data$outlier <- pca_data$distance > threshold

# PCA Scatter Plot
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = outlier)) +
  geom_jitter(width = 0.02, height = 0.02, size = 1.8, alpha = 0.6) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  labs(x = "PC1", y = "PC2", color = "Outlier") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank()) +
  ggrepel::geom_text_repel(
    data = subset(pca_data, outlier),
    aes(label = Individual),
    size = 3, color = "black"
  )

# Save PCA Plot
ggsave("results/PCA_plot_refined.png", plot = pca_plot, width = 6, height = 5, dpi = 300)
cat("📊 PCA plot saved to results/PCA_plot_refined.png\n")

# -----------------------------
# STEP 3: VARIANCE EXPLAINED
# -----------------------------
eigenvalues <- scan("pf_snp_pca.eigenval")
variance_explained <- (eigenvalues / sum(eigenvalues)) * 100
cat("PC1 explains:", round(variance_explained[1], 2), "%\n")
cat("PC2 explains:", round(variance_explained[2], 2), "%\n")

# -----------------------------
# STEP 4: HIERARCHICAL CLUSTERING
# -----------------------------
distance_matrix <- dist(pca_data[, 3:12], method = "euclidean")
hc <- hclust(distance_matrix, method = "ward.D2")

# Save dendrogram (base R)
png("results/Hierarchical_Clustering_Dendrogram.png", width = 8, height = 6, units = "in", res = 300)
plot(hc, labels = pca_data$Individual, main = "", cex = 0.3)
dev.off()
cat("🌳 Dendrogram saved to results/Hierarchical_Clustering_Dendrogram.png\n")

# -----------------------------
# STEP 5: GGPLOT2 DENDROGRAM
# -----------------------------
dendro_data <- as.dendrogram(hc)
dendro_df <- ggdendro::dendro_data(dendro_data)

dendrogram_plot <- ggplot() +
  geom_segment(data = dendro_df$segments,
               aes(x = x, y = y, xend = xend, yend = yend),
               color = "blue", size = 0.5) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank()) +
  labs(title = "Hierarchical Clustering Dendrogram")

ggsave("results/Dendrogram_Plot.png", plot = dendrogram_plot, width = 8, height = 2, dpi = 300)

# -----------------------------
# STEP 6: STATISTICAL TESTS
# -----------------------------
pca_data$outlier_group <- ifelse(pca_data$outlier, "Outlier", "Main Cluster")
anova_result <- aov(distance ~ outlier_group, data = pca_data)
cat("\n📊 ANOVA Result:\n")
print(summary(anova_result))

# -----------------------------
# STEP 7: COPHENETIC CORRELATION
# -----------------------------
cophenetic_dist <- cophenetic(hc)
ccc <- cor(distance_matrix, cophenetic_dist)
cat("\n🔢 Cophenetic Correlation Coefficient:", round(ccc, 4), "\n")
