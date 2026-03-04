#!/bin/bash
# ===============================================================
# Script: 04_filter_and_summarize_snps.sh
# Purpose: Filter SNPs for resistance genes of interest and
#          generate summary statistics including counts per gene,
#          per year, and amino acid change frequencies.
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
#   - filtered_snps_resistance_genes.tsv     SNPs for genes of interest
#   - snp_count_per_gene.txt                 SNP count per gene
#   - snp_count_per_year_per_gene.tsv        SNP count per year per gene
#   - amino_acid_change_freq.tsv             Amino acid change frequencies
#
# GENES OF INTEREST:
#   PF3D7_0709000  pfcrt   (chloroquine resistance transporter)
#   PF3D7_0523000  pfmdr1  (multidrug resistance 1)
#   PF3D7_1408000  pfdhps  (dihydropteroate synthase)
#   PF3D7_1408100  pfdhps  (dihydropteroate synthase, adjacent)
#   PF3D7_0417200  pfdhfr  (dihydrofolate reductase)
#   PF3D7_0810800  pfpm2   (plasmepsin 2, copy number variation)
#   PF3D7_1343700  pfk13   (kelch 13, artemisinin resistance)
#
# USAGE:
#   bash 04_filter_and_summarize_snps.sh
# ===============================================================

set -euo pipefail

# ---------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------
INPUT="PLINK_SNPs_summary_per_year_fixed.tsv"

# Output files
FILTERED="filtered_snps_resistance_genes.tsv"
SNP_COUNT="snp_count_per_gene.txt"
YEARLY_COUNT="snp_count_per_year_per_gene.tsv"
AA_FREQ="amino_acid_change_freq.tsv"

# Genes of interest (column 7 = Gene_ID in input file)
GENES=(
    "PF3D7_0709000"   # pfcrt
    "PF3D7_0523000"   # pfmdr1
    "PF3D7_1408000"   # pfdhps
    "PF3D7_1408100"   # pfdhps (adjacent)
    "PF3D7_0417200"   # pfdhfr
    "PF3D7_0810800"   # pfpm2
    "PF3D7_1343700"   # pfk13
)

# ---------------------------------------------------------------
# VALIDATE INPUT
# ---------------------------------------------------------------
if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: Input file not found: $INPUT"
    echo "       Please ensure PLINK_SNPs_summary_per_year_fixed.tsv"
    echo "       is in the current directory before running this script."
    exit 1
fi

echo "Input file: $INPUT"
echo "Total rows (excl. header): $(tail -n +2 "$INPUT" | wc -l)"

# ---------------------------------------------------------------
# STEP 1: FILTER SNPs BY GENE OF INTEREST
# ---------------------------------------------------------------
echo ""
echo ">>> Filtering SNPs for $(echo ${#GENES[@]}) resistance genes..."

# Build awk condition dynamically from GENES array
AWK_COND=$(printf '($7 == "%s") || ' "${GENES[@]}")
AWK_COND=${AWK_COND::-4}  # remove trailing " || "

# Write header + filtered rows
head -n 1 "$INPUT" > "$FILTERED"
awk -F"\t" "$AWK_COND { print }" "$INPUT" >> "$FILTERED"

FILTERED_COUNT=$(tail -n +2 "$FILTERED" | wc -l)
echo "    Retained $FILTERED_COUNT SNP records across ${#GENES[@]} genes"
echo "    Saved to: $FILTERED"

# ---------------------------------------------------------------
# STEP 2: SNP COUNT PER GENE
# ---------------------------------------------------------------
echo ""
echo ">>> Counting SNPs per gene..."

{
    echo -e "Gene\tSNP_Count"
    cut -f7 "$FILTERED" | tail -n +2 | sort | uniq -c | sort -nr | \
    awk '{print $2 "\t" $1}'
} > "$SNP_COUNT"

echo "    Results:"
column -t "$SNP_COUNT"
echo "    Saved to: $SNP_COUNT"

# ---------------------------------------------------------------
# STEP 3: SNP COUNT PER YEAR PER GENE
# ---------------------------------------------------------------
echo ""
echo ">>> Counting SNPs per year per gene..."

{
    echo -e "Year\tGene\tSNP_Count"
    tail -n +2 "$FILTERED" | cut -f1,7 | sort | uniq -c | \
    awk '{print $2 "\t" $3 "\t" $1}'
} > "$YEARLY_COUNT"

echo "    Results:"
column -t "$YEARLY_COUNT"
echo "    Saved to: $YEARLY_COUNT"

# ---------------------------------------------------------------
# STEP 4: AMINO ACID CHANGE FREQUENCY PER GENE
# ---------------------------------------------------------------
echo ""
echo ">>> Calculating amino acid change frequencies..."

{
    echo -e "Gene\tAmino_Acid_Change\tCount"
    tail -n +2 "$FILTERED" | cut -f7,8 | sort | uniq -c | \
    sort -k2,2 -k1,1nr | \
    awk '{print $2 "\t" $3 "\t" $1}'
} > "$AA_FREQ"

echo "    Results:"
column -t "$AA_FREQ"
echo "    Saved to: $AA_FREQ"

# ---------------------------------------------------------------
# SUMMARY
# ---------------------------------------------------------------
echo ""
echo ">>> All done. Output files:"
echo "    $FILTERED      — filtered SNP records"
echo "    $SNP_COUNT     — SNP count per gene"
echo "    $YEARLY_COUNT  — SNP count per year per gene"
echo "    $AA_FREQ       — amino acid change frequencies"
