# ===============================================================
# Script:  03_snp_trend_analysis.R
# Purpose: SNP allele frequency trend analysis for key
#          antimalarial resistance markers in P. falciparum
#          Upper East Region, Ghana (2009-2018)
# Study:   Longitudinal Genomic Surveillance of Antimalarial
#          Drug Resistance-Associated SNPs, Upper East Ghana
# Author:  Hayford Osei Offei
# Date:    2025
# ===============================================================
#
# INPUT FILES REQUIRED (place in working directory):
#   - Pf7_drug_resistance_marker_genotypes.txt  (MalariaGEN Pf7)
#   - Pf7_samples.txt                           (MalariaGEN Pf7)
#   - filtered_samples_YYYY.csv                 (one per year, 2009-2018)
#   - combined_haplotypes_NposAA_summary.csv    (for haplotype analysis)
#
# OUTPUT FILES (saved to results/):
#   MANUSCRIPT FIGURES:
#   - crt_CVIET_CVMNK_trend_plot_fixed.png      (Figure: pfcrt haplotypes)
#   - double_triple_mutant_plot.png             (Figure: pfdhfr haplotypes)
#   - SNP_trend_kelch13_*.png                   (Figure: pfk13 trends)
#   - SNP_trend_mdr1_*.png                      (Figure: pfmdr1 CNV trends)
#
#   SUPPLEMENTARY (all other SNPs):
#   - SNP_trend_[snp_name].png                  (one per SNP marker)
#
#   DATA OUTPUTS:
#   - SNP_trend_pvalues.csv                     (chi-square trend p-values)
#   - crt_CVIET_CVMNK_trend_data.csv
#   - double_triple_mutant_data.csv
#
# DEPENDENCIES:
#   - R packages: dplyr, tidyr, tidyverse, ggplot2
# ===============================================================


# ---------------------------------------------------------------
# SETUP
# ---------------------------------------------------------------
setwd("E:/Malaria_Project/New documents")  # Update path as needed

dir.create("results", showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tidyverse)
})


# ---------------------------------------------------------------
# STEP 1: LOAD AND MERGE HAPLOTYPE METADATA WITH YEARLY DATA
# ---------------------------------------------------------------
# Load MalariaGEN Pf7 drug resistance marker genotypes
haplo_df <- read.delim(
  "Pf7_drug_resistance_marker_genotypes.txt",
  stringsAsFactors = FALSE
)
names(haplo_df)[1] <- "sample_id"

cat("Loaded haplotype metadata:", nrow(haplo_df), "samples\n")

# Merge with filtered per-year sample files
years <- 2009:2018

for (year in years) {
  input_file  <- sprintf("filtered_samples_%d.csv", year)
  output_file <- sprintf("merged_samples_%d.csv", year)

  if (!file.exists(input_file)) {
    cat("Skipping", input_file, "- file not found\n")
    next
  }

  read.csv(input_file, stringsAsFactors = FALSE) %>%
    inner_join(haplo_df, by = "sample_id") %>%
    write.csv(output_file, row.names = FALSE)

  cat("Merged and saved:", output_file, "\n")
}


# ---------------------------------------------------------------
# STEP 2: LOAD AND COMBINE ALL MERGED YEARLY DATA
# ---------------------------------------------------------------
merged_list <- lapply(years, function(y) {
  file <- sprintf("merged_samples_%d.csv", y)
  if (file.exists(file)) {
    read.csv(file, stringsAsFactors = FALSE) %>%
      mutate(year = y)
  } else NULL
})

all_data <- bind_rows(merged_list)
if (nrow(all_data) == 0) stop("No merged data found. Please check input files.")
cat("Combined dataset:", nrow(all_data), "rows across", n_distinct(all_data$year), "years\n")


# ---------------------------------------------------------------
# STEP 3: ALL-SNP TREND PLOTS (supplementary figures)
# Generates one bar chart per SNP marker showing yearly
# allele frequency trends. Covers pfcrt, pfdhfr, pfdhps,
# pfmdr1, pfk13, and pfpm2 loci.
# ---------------------------------------------------------------
snp_cols <- grep(
  "^(crt_|dhfr_|dhps_|kelch13_|mdr1_|pm2_)",
  names(all_data),
  value = TRUE
)
all_data[snp_cols] <- lapply(all_data[snp_cols], as.character)

snp_data_long <- all_data %>%
  select(year, all_of(snp_cols)) %>%
  pivot_longer(
    cols      = -year,
    names_to  = "snp",
    values_to = "allele"
  ) %>%
  filter(allele != "" & allele != "." & !is.na(allele)) %>%
  mutate(allele = ifelse(allele == "-", "Untyped", allele))

# Summarise allele counts and percentages per year
year_totals <- all_data %>%
  group_by(year) %>%
  summarise(total_samples = n(), .groups = "drop")

snp_summary <- snp_data_long %>%
  group_by(year, snp, allele) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(year_totals, by = "year") %>%
  mutate(percent = (count / total_samples) * 100)

# Generate one plot per SNP
for (snp_name in unique(snp_summary$snp)) {
  snp_plot_data <- snp_summary %>% filter(snp == snp_name)

  p <- ggplot(snp_plot_data,
              aes(x = factor(year), y = percent, fill = allele)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(
      aes(label = sprintf("%.1f%%", percent)),
      position  = position_stack(vjust = 0.5),
      color     = "white",
      fontface  = "bold",
      size      = 3.5
    ) +
    scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
    labs(
      x    = "Year",
      y    = "Percentage (%)",
      fill = "Allele"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x       = element_text(angle = 45, hjust = 1),
      legend.position   = "right",
      legend.title      = element_text(face = "bold"),
      panel.grid.major  = element_line(color = "gray85"),
      panel.grid.minor  = element_blank()
    )

  plot_file <- paste0("results/SNP_trend_", snp_name, ".png")
  ggsave(plot_file, plot = p, width = 8, height = 5, dpi = 300)
  cat("Saved:", plot_file, "\n")
}


# ---------------------------------------------------------------
# STEP 4: MANUSCRIPT FIGURE — pfcrt CVIET vs CVMNK HAPLOTYPES
# Plots yearly frequency of chloroquine resistance (CVIET) vs
# wild-type (CVMNK) pfcrt haplotypes (2009-2018).
# ---------------------------------------------------------------
crt_cols <- grep("^crt_", names(all_data), value = TRUE)
if (length(crt_cols) == 0) stop("No crt columns found in dataset.")

crt_df <- all_data %>%
  select(year, sample_id, all_of(crt_cols)) %>%
  pivot_longer(
    cols      = all_of(crt_cols),
    names_to  = "crt_marker",
    values_to = "haplotype"
  ) %>%
  filter(!is.na(haplotype) & haplotype != "" & haplotype != ".") %>%
  filter(grepl("CVIET|CVMNK", haplotype, ignore.case = TRUE)) %>%
  mutate(
    haplotype_group = case_when(
      grepl("CVIET", haplotype, ignore.case = TRUE) ~ "CVIET (Resistant)",
      grepl("CVMNK", haplotype, ignore.case = TRUE) ~ "CVMNK (Wild-type)",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(haplotype_group))

crt_year_totals <- crt_df %>%
  group_by(year) %>%
  summarise(total_samples = n(), .groups = "drop")

crt_summary <- crt_df %>%
  group_by(year, haplotype_group) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(crt_year_totals, by = "year") %>%
  mutate(actual_percent = (count / total_samples) * 100)

crt_colors <- c(
  "CVIET (Resistant)" = "#E67E22",  # orange
  "CVMNK (Wild-type)" = "#3498DB"   # blue
)

if (nrow(crt_summary) > 0) {
  p_crt <- ggplot(
    crt_summary,
    aes(x = factor(year), y = actual_percent, fill = haplotype_group)
  ) +
    geom_col(
      position = position_dodge(width = 0.8),
      color    = "black",
      width    = 0.7,
      alpha    = 0.9
    ) +
    geom_text(
      aes(label = paste0(round(actual_percent, 1), "%")),
      position = position_dodge(width = 0.8),
      vjust    = -0.5,
      hjust    = 0.5,
      size     = 4.0,
      fontface = "bold",
      color    = "#1B2631"
    ) +
    scale_fill_manual(values = crt_colors, name = "pfcrt Haplotype") +
    scale_y_continuous(limits = c(0, 115), expand = c(0, 0)) +
    labs(
      x = "Year",
      y = "Percentage (%)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.title        = element_text(face = "bold", size = 13),
      axis.text         = element_text(size = 12, color = "#2C3E50"),
      axis.text.x       = element_text(angle = 45, hjust = 1),
      panel.grid.minor  = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position   = "top",
      legend.title      = element_text(face = "bold"),
      legend.text       = element_text(size = 11)
    )

  ggsave(
    "results/crt_CVIET_CVMNK_trend_plot_fixed.png",
    plot   = p_crt,
    width  = 9,
    height = 6,
    dpi    = 300
  )
  write.csv(crt_summary, "results/crt_CVIET_CVMNK_trend_data.csv",
            row.names = FALSE)
  cat("pfcrt haplotype plot saved\n")
} else {
  cat("No pfcrt CVIET/CVMNK data found\n")
}


# ---------------------------------------------------------------
# STEP 5: MANUSCRIPT FIGURE — pfdhfr DOUBLE & TRIPLE MUTANT
#         HAPLOTYPES
# Classifies and plots yearly frequency of double mutants
# (C59R+S108N, N51I+S108N) and triple mutant (N51I+C59R+S108N)
# antifolate resistance haplotypes.
# ---------------------------------------------------------------
df_haplo <- read.csv(
  "combined_haplotypes_NposAA_summary.csv",
  stringsAsFactors = FALSE
)

# Helper: check if a string contains all specified patterns
contains_all <- function(str, patterns) {
  all(sapply(patterns, function(p) grepl(p, str, ignore.case = TRUE)))
}

# Classify haplotype groups (order-insensitive matching)
mutant_df <- df_haplo %>%
  rowwise() %>%
  mutate(
    haplotype_group = case_when(
      contains_all(combined_haplotype, c("C59R", "S108N")) &
        !grepl("N51I", combined_haplotype, ignore.case = TRUE) ~
        "C59R+S108N (Double)",
      contains_all(combined_haplotype, c("N51I", "S108N")) &
        !grepl("C59R", combined_haplotype, ignore.case = TRUE) ~
        "N51I+S108N (Double)",
      contains_all(combined_haplotype, c("N51I", "C59R", "S108N")) ~
        "N51I+C59R+S108N (Triple)",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  filter(!is.na(haplotype_group))

# Calculate percentages
haplo_year_totals <- df_haplo %>%
  group_by(year) %>%
  summarise(total_samples = sum(count), .groups = "drop")

mutant_df <- mutant_df %>%
  left_join(haplo_year_totals, by = "year") %>%
  mutate(actual_percent = (count / total_samples) * 100) %>%
  group_by(year, haplotype_group) %>%
  summarise(actual_percent = sum(actual_percent, na.rm = TRUE),
            .groups = "drop")

haplo_colors <- c(
  "C59R+S108N (Double)"      = "#2ECC71",  # green
  "N51I+S108N (Double)"      = "#E67E22",  # orange
  "N51I+C59R+S108N (Triple)" = "#3498DB"   # blue
)

if (nrow(mutant_df) > 0) {
  p_haplo <- ggplot(
    mutant_df,
    aes(x = factor(year), y = actual_percent, fill = haplotype_group)
  ) +
    geom_col(
      position = position_dodge(width = 0.8),
      color    = "black",
      width    = 0.7,
      alpha    = 0.9
    ) +
    geom_text(
      aes(label = paste0(round(actual_percent, 1), "%")),
      position = position_dodge(width = 0.8),
      vjust    = -0.3,
      hjust    = 0.5,
      size     = 3.5,
      fontface = "bold",
      color    = "#1B2631"
    ) +
    scale_fill_manual(values = haplo_colors, name = "Haplotype Group") +
    scale_y_continuous(limits = c(0, 115), expand = c(0, 0)) +
    labs(
      x = "Year",
      y = "Percentage (%)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.title         = element_text(face = "bold", size = 13),
      axis.text          = element_text(size = 12, color = "#2C3E50"),
      axis.text.x        = element_text(angle = 45, hjust = 1),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position    = "top",
      legend.title       = element_text(face = "bold"),
      legend.text        = element_text(size = 11)
    )

  ggsave(
    "results/double_triple_mutant_plot.png",
    plot   = p_haplo,
    width  = 9,
    height = 6,
    dpi    = 300
  )
  write.csv(mutant_df, "results/double_triple_mutant_data.csv",
            row.names = FALSE)
  cat("pfdhfr haplotype plot saved\n")
} else {
  cat("No matching haplotypes found in dataset\n")
}


# ---------------------------------------------------------------
# STEP 6: TREND TESTS (Pearson's chi-square)
# Tests whether allele frequencies changed significantly
# over time for each SNP marker. p < 0.05 = significant trend.
# ---------------------------------------------------------------
trend_results <- snp_data_long %>%
  group_by(snp) %>%
  group_modify(~ {
    tbl <- table(.x$year, .x$allele)
    if (nrow(tbl) > 1 && ncol(tbl) > 1) {
      allele_counts <- tbl[, 1]
      total_counts  <- rowSums(tbl)
      if (length(allele_counts) == length(total_counts) &&
          all(total_counts > 0)) {
        p_val <- prop.trend.test(allele_counts, total_counts)$p.value
      } else {
        p_val <- NA
      }
    } else {
      p_val <- NA
    }
    tibble(p_value = p_val)
  }) %>%
  ungroup() %>%
  mutate(
    significant = case_when(
      is.na(p_value)  ~ "Insufficient data",
      p_value < 0.05  ~ "Yes",
      TRUE            ~ "No"
    )
  ) %>%
  arrange(p_value)

write.csv(trend_results, "results/SNP_trend_pvalues.csv", row.names = FALSE)
cat("SNP trend p-values saved to results/SNP_trend_pvalues.csv\n")
cat("\nSignificant SNP trends (p < 0.05):\n")
print(trend_results %>% filter(significant == "Yes"))
