#!/usr/bin/env Rscript

# ===============================================================
# Script: snp_analysis_pipeline.R
# Purpose: End-to-end pipeline for merging haplotype metadata,
#          combining yearly data, summarizing SNP allele frequencies,
#          visualizing trends, and performing p-value trend tests.
# Author: [Your Name]
# Usage: Rscript snp_analysis_pipeline.R
# ===============================================================

library(dplyr)
library(tidyverse)

# --- 1. Merge haplotype metadata with filtered yearly data ---
haplo_df <- read.delim("Pf7_drug_resistance_marker_genotypes.txt", stringsAsFactors = FALSE)
names(haplo_df)[1] <- "sample_id"

years <- 2009:2018
for (year in years) {
  input_file  <- sprintf("filtered_samples_%d.csv", year)
  output_file <- sprintf("merged_samples_%d.csv", year)

  if (!file.exists(input_file)) {
    cat("⚠ Skipping", input_file, "- file not found\n")
    next
  }

  read.csv(input_file, stringsAsFactors = FALSE) %>%
    inner_join(haplo_df, by = "sample_id") %>%
    write.csv(output_file, row.names = FALSE)

  cat("✔ Merged and saved:", output_file, "\n")
}

# --- 2. Load and combine merged yearly data ---
merged_list <- lapply(years, function(y) {
  file <- sprintf("merged_samples_%d.csv", y)
  if (file.exists(file)) {
    read.csv(file, stringsAsFactors = FALSE) %>%
      mutate(year = y)
  } else NULL
})

all_data <- bind_rows(merged_list)
if (nrow(all_data) == 0) stop("No merged data found. Please check input files.")

# --- 3. Identify SNP columns & reshape ---
snp_cols <- grep("^(crt_|dhfr_|dhps_|kelch13_|mdr1_|pm2_)", names(all_data), value = TRUE)
all_data[snp_cols] <- lapply(all_data[snp_cols], as.character)

snp_data_long <- all_data %>%
  select(year, all_of(snp_cols)) %>%
  pivot_longer(cols = -year, names_to = "snp", values_to = "allele") %>%
  filter(allele != "" & allele != "." & !is.na(allele))

# --- 4. Summarize allele counts and percentages ---
year_totals <- all_data %>%
  group_by(year) %>%
  summarise(total_samples = n(), .groups = "drop")

snp_summary <- snp_data_long %>%
  group_by(year, snp, allele) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(year_totals, by = "year") %>%
  mutate(percent = (count / total_samples) * 100)

# --- 5. Plot allele frequency trends ---
plot_file <- "SNP_allele_trends.png"
png(plot_file, width = 1600, height = 900)
ggplot(snp_summary, aes(x = factor(year), y = percent, fill = allele)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ snp, scales = "free_y", ncol = 3) +
  ylab("Allele Frequency (%)") + xlab("Year") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
dev.off()
cat("✔ Allele frequency plot saved to", plot_file, "\n")

# --- 6. Perform trend tests (p-values) ---
trend_results <- snp_data_long %>%
  group_by(snp) %>%
  group_modify(~{
    tbl <- table(.x$year, .x$allele)
    if (nrow(tbl) > 1 && ncol(tbl) > 1) {
      allele_counts <- tbl[, 1]              # first allele only
      total_counts <- rowSums(tbl)
      if (length(allele_counts) == length(total_counts) && all(total_counts > 0)) {
        p_val <- prop.trend.test(allele_counts, total_counts)$p.value
      } else p_val <- NA
    } else p_val <- NA
    tibble(p_value = p_val)
  }) %>%
  ungroup() %>%
  mutate(significant = case_when(
    is.na(p_value) ~ "Not enough data",
    p_value < 0.05 ~ "Yes",
    TRUE ~ "No"
  ))

write.csv(trend_results, "SNP_trend_pvalues.csv", row.names = FALSE)
cat("✔ SNP trend p-values saved to SNP_trend_pvalues.csv\n")
