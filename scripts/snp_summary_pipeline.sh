#!/bin/bash
# ===============================================================
# Script: filter_and_summarize_snps.sh
# Purpose: Filter SNPs for genes of interest and generate summary 
#          statistics including counts per gene, per year, and 
#          amino acid change frequencies.
# Author: [Hayford Osei Offei]
# Date: [18-09-2025]
# ===============================================================

# -----------------------------
# INPUT / OUTPUT FILES
# -----------------------------
INPUT="PLINK_SNPs_summary_per_year_fixed.tsv"
FILTERED="filtered_snps_MYGENES.tsv"
SNP_COUNT="snp_count_per_gene_MYGENES.txt"
YEARLY_COUNT="snp_count_per_year_per_gene_MYGENES.tsv"
AA_FREQ="amino_acid_change_freq_MYGENES.tsv"

# -----------------------------
# GENES OF INTEREST
# -----------------------------
GENES=(
    "PF3D7_0709000"
    "PF3D7_0523000"
    "PF3D7_1408000"
    "PF3D7_1408100"
    "PF3D7_0417200"
    "PF3D7_0810800"
    "PF3D7_1343700"
)

# -----------------------------
# BUILD AWK CONDITION
# -----------------------------
AWK_COND=$(printf '($7 == "%s") || ' "${GENES[@]}")
AWK_COND=${AWK_COND::-4}  # remove trailing " || "

# -----------------------------
# STEP 1: FILTER SNPs BY GENE
# -----------------------------
echo "🔎 Filtering SNPs for genes of interest..."
head -n 1 "$INPUT" > "$FILTERED"  # keep header
awk -F"\t" "$AWK_COND { print }" "$INPUT" >> "$FILTERED"
echo "✔ Filtered SNPs saved to $FILTERED"

# -----------------------------
# STEP 2: SNP COUNT PER GENE
# -----------------------------
echo -e "\n📊 SNP count per gene:"
{
    echo -e "Gene\tSNP_Count"
    cut -f7 "$FILTERED" | tail -n +2 | sort | uniq -c | sort -nr | \
    awk '{print $2 "\t" $1}'
} > "$SNP_COUNT"
cat "$SNP_COUNT"

# -----------------------------
# STEP 3: SNP COUNT PER YEAR PER GENE
# -----------------------------
echo -e "\n📊 SNP count per year per gene:"
{
    echo -e "Year\tGene\tSNP_Count"
    tail -n +2 "$FILTERED" | cut -f1,7 | sort | uniq -c | \
    awk '{print $2 "\t" $3 "\t" $1}'
} > "$YEARLY_COUNT"
column -t "$YEARLY_COUNT"

# -----------------------------
# STEP 4: AMINO ACID CHANGE FREQUENCY
# -----------------------------
echo -e "\n📊 Amino acid change frequency per gene:"
{
    echo -e "Gene\tAmino_Acid_Change\tCount"
    tail -n +2 "$FILTERED" | cut -f7,8 | sort | uniq -c | \
    sort -k2,2 -k1,1nr | \
    awk '{print $2 "\t" $3 "\t" $1}'
} > "$AA_FREQ"
column -t "$AA_FREQ"

# -----------------------------
# SUMMARY MESSAGE
# -----------------------------
echo -e "\n✅ Done. Output files created:"
echo " - $FILTERED"
echo " - $SNP_COUNT"
echo " - $YEARLY_COUNT"
echo " - $AA_FREQ"


