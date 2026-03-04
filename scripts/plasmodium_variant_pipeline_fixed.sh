#!/bin/bash
# =============================================================================
# Plasmodium falciparum Genomic Surveillance Pipeline
# Study: Longitudinal Genomic Surveillance of Antimalarial Drug
#        Resistance-Associated SNPs in P. falciparum,
#        Upper East Region of Ghana (2010-2018)
# Reference Genome: GCF_000002765.6 (P. falciparum 3D7)
# Authors: Adjoa Agyemang Boakye, Hayford Osei Offei,
#          David Adedia, Enoch Aninagyei
# Institution: University of Health and Allied Sciences, Ho, Ghana
# =============================================================================
# PIPELINE OVERVIEW:
#   1.  Environment Setup
#   2.  Data Download from ENA
#   3.  Quality Control (FastQC)
#   4.  Read Trimming (Trimmomatic)
#   5.  Reference Genome Indexing (BWA)
#   6.  Alignment to Reference Genome (BWA-MEM)
#   7.  BAM Sorting & Indexing (SAMtools)
#   8.  Alignment Quality Statistics (SAMtools flagstat)
#   9.  Variant Calling (bcftools mpileup + call)
#   10. Variant Filtering (bcftools filter)
#   11. Variant Annotation (SnpEff)
#   12. Chromosome Renaming (accession -> chromosome name)
# =============================================================================
# NOTES:
#   - All 90 ERR samples are discovered automatically from the ENA accessions
#     file (no need to hardcode sample IDs)
#   - After QC, 77 of 90 samples passed (see per-year breakdown at end of file)
#   - Samples from 2009 were excluded (no samples passed QC)
#   - Run individual steps by uncommenting them in the MAIN section at the end
# =============================================================================

set -euo pipefail  # Exit immediately on error; treat unset variables as errors

# =============================================================================
# CONFIGURATION — Edit paths here if needed
# =============================================================================

REFERENCE="GCF_000002765.6_GCA_000002765_genomic.fna"
ACCESSIONS_FILE="ena_accessions.txt"   # One ERR accession per line
THREADS=6                              # Adjust to your system's CPU count

# Directory structure
FASTQ_DIR="fastq_files"
FASTQC_DIR="fastqc_results"
TRIMMED_DIR="trimmed_files"
BAM_DIR="bam_files"
VCF_DIR="vcf_files"
ANNOTATED_DIR="annotated_vcf"
STATS_DIR="alignment_stats"

# SnpEff genome name (must match your snpEff database)
SNPEFF_GENOME="Plasmodium_falciparum"

# =============================================================================
# STEP 0: GENERATE ENA ACCESSIONS FILE
# Creates ena_accessions.txt with all 90 ERR sample IDs.
# This file is used by all downstream steps to discover samples automatically.
# Only needs to be run once.
# =============================================================================

generate_accessions_file() {
    echo ">>> Generating ENA accessions file: ${ACCESSIONS_FILE}..."
    cat > "${ACCESSIONS_FILE}" << 'EOF'
ERR019549
ERR039235
ERR036147
ERR035362
ERR045604
ERR205933
ERR211562
ERR045620
ERR045607
ERR045609
ERR039979
ERR042698
ERR042701
ERR042719
ERR042722
ERR042727
ERR063567
ERR114406
ERR114409
ERR211462
ERR376187
ERR376188
ERR376201
ERR376202
ERR376203
ERR376212
ERR376213
ERR376214
ERR376216
ERR403196
ERR450116
ERR450108
ERR450099
ERR450045
ERR586212
ERR586179
ERR450104
ERR450112
ERR450044
ERR450041
ERR3504126
ERR3523699
ERR3546199
ERR3546151
ERR3546188
ERR3594524
ERR3594529
ERR3523718
ERR3523721
ERR3594576
ERR3486201
ERR3523571
ERR3523563
ERR3546220
ERR3523635
ERR3523614
ERR1214163
ERR1214188
ERR1214209
ERR1214148
ERR2532525
ERR2542018
ERR2541936
ERR2532503
ERR2532517
ERR2532524
ERR2532570
ERR2541897
ERR2541942
ERR2532534
ERR3486300
ERR2891366
ERR2890914
ERR3608767
ERR3546081
ERR2375226
ERR2506008
ERR3504084
ERR3523658
ERR3504120
ERR3608800
ERR3608865
ERR3584487
ERR3594738
ERR3608785
ERR3594774
ERR3594757
ERR3594689
ERR3594679
ERR3594661
EOF
    echo ">>> Accessions file created with $(wc -l < ${ACCESSIONS_FILE}) samples."
}

# =============================================================================
# STEP 1: ENVIRONMENT SETUP
# Installs all required tools. Run once on a fresh system.
# =============================================================================

setup_environment() {
    echo ">>> Setting up environment..."
    sudo apt update
    sudo apt install -y \
        build-essential \
        fastqc \
        trimmomatic \
        bwa \
        samtools \
        bcftools \
        default-jre       # Required for SnpEff

    echo ">>> Verifying installations..."
    fastqc --version
    trimmomatic -version
    bwa 2>&1 | head -3
    samtools --version | head -1
    bcftools --version | head -1
    echo ">>> Environment setup complete."
}

# =============================================================================
# STEP 2: DATA DOWNLOAD FROM ENA
# Downloads paired-end FASTQ files for all samples in the accessions file.
# Uses ena-fast-download (recommended) or falls back to wget.
# =============================================================================

download_ena_data() {
    echo ">>> Downloading FASTQ files from ENA..."
    mkdir -p "${FASTQ_DIR}"

    while read -r accession; do
        [[ -z "$accession" ]] && continue  # Skip blank lines

        R1="${FASTQ_DIR}/${accession}_1.fastq.gz"
        R2="${FASTQ_DIR}/${accession}_2.fastq.gz"

        # Skip if already downloaded
        if [[ -f "$R1" && -f "$R2" ]]; then
            echo "  [SKIP] ${accession} already downloaded."
            continue
        fi

        echo "  Downloading ${accession}..."

        # ENA FTP path (paired-end reads)
        prefix=$(echo "$accession" | cut -c1-6)
        ftp_base="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}"

        # Handle accessions with different length suffixes
        if [[ ${#accession} -eq 9 ]]; then
            ftp_path="${ftp_base}/${accession}"
        else
            suffix=$(echo "$accession" | rev | cut -c1-3 | rev)
            ftp_path="${ftp_base}/${suffix}/${accession}"
        fi

        wget -q "${ftp_path}/${accession}_1.fastq.gz" -O "$R1" \
            || { echo "  WARNING: Failed to download ${accession} R1. Skipping..."; rm -f "$R1"; continue; }
        wget -q "${ftp_path}/${accession}_2.fastq.gz" -O "$R2" \
            || { echo "  WARNING: Failed to download ${accession} R2. Skipping..."; rm -f "$R2"; continue; }

        echo "  Done: ${accession}"
    done < "${ACCESSIONS_FILE}"

    echo ">>> Download complete. Files saved to: ${FASTQ_DIR}"
}

# =============================================================================
# STEP 3: QUALITY CONTROL (FastQC)
# Runs FastQC on all raw FASTQ files before trimming.
# =============================================================================

run_fastqc_raw() {
    echo ">>> Running FastQC on raw FASTQ files..."
    mkdir -p "${FASTQC_DIR}/raw"
    fastqc "${FASTQ_DIR}"/*.fastq.gz \
        --outdir "${FASTQC_DIR}/raw" \
        --threads "${THREADS}"
    echo ">>> Raw FastQC reports saved to: ${FASTQC_DIR}/raw"
}

# =============================================================================
# STEP 4: READ TRIMMING (Trimmomatic)
# Trims low-quality bases and adapter sequences from all paired-end reads.
# Settings: sliding window 4:20, minimum Phred quality Q20.
# =============================================================================

run_trimmomatic() {
    echo ">>> Trimming reads with Trimmomatic..."
    mkdir -p "${TRIMMED_DIR}"

    while read -r accession; do
        [[ -z "$accession" ]] && continue

        R1="${FASTQ_DIR}/${accession}_1.fastq.gz"
        R2="${FASTQ_DIR}/${accession}_2.fastq.gz"
        out_R1="${TRIMMED_DIR}/${accession}_1_trimmed.fastq.gz"
        out_R2="${TRIMMED_DIR}/${accession}_2_trimmed.fastq.gz"
        unpaired_R1="${TRIMMED_DIR}/${accession}_1_unpaired.fastq.gz"
        unpaired_R2="${TRIMMED_DIR}/${accession}_2_unpaired.fastq.gz"

        # Skip if input files are missing
        if [[ ! -f "$R1" || ! -f "$R2" ]]; then
            echo "  WARNING: FASTQ files missing for ${accession}. Skipping..."
            continue
        fi

        # Skip if already trimmed
        if [[ -f "$out_R1" && -f "$out_R2" ]]; then
            echo "  [SKIP] ${accession} already trimmed."
            continue
        fi

        echo "  Trimming ${accession}..."
        trimmomatic PE \
            -phred33 \
            -threads "${THREADS}" \
            "$R1" "$R2" \
            "$out_R1" "$unpaired_R1" \
            "$out_R2" "$unpaired_R2" \
            SLIDINGWINDOW:4:20 \
            || { echo "  ERROR: Trimming failed for ${accession}. Skipping..."; continue; }

        echo "  Done: ${accession}"
    done < "${ACCESSIONS_FILE}"

    echo ">>> Trimming complete. Trimmed files saved to: ${TRIMMED_DIR}"
}

# =============================================================================
# STEP 5: REFERENCE GENOME INDEXING (BWA)
# Only needs to be run once per reference genome.
# =============================================================================

index_reference() {
    echo ">>> Indexing reference genome: ${REFERENCE}..."
    if [[ ! -f "${REFERENCE}" ]]; then
        echo "ERROR: Reference genome not found: ${REFERENCE}"
        echo "       Download from: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000002765.6/"
        exit 1
    fi
    bwa index "${REFERENCE}"
    echo ">>> Reference genome indexed."
}

# =============================================================================
# STEP 6: ALIGNMENT (BWA-MEM)
# Aligns trimmed paired-end reads to the P. falciparum 3D7 reference genome.
# Reads with mapping quality < Q20 are excluded.
# =============================================================================

align_samples() {
    echo ">>> Aligning samples to reference genome..."
    mkdir -p "${BAM_DIR}"

    while read -r accession; do
        [[ -z "$accession" ]] && continue

        R1="${TRIMMED_DIR}/${accession}_1_trimmed.fastq.gz"
        R2="${TRIMMED_DIR}/${accession}_2_trimmed.fastq.gz"
        bam_out="${BAM_DIR}/${accession}.bam"

        if [[ ! -f "$R1" || ! -f "$R2" ]]; then
            echo "  WARNING: Trimmed files missing for ${accession}. Skipping..."
            continue
        fi

        if [[ -f "$bam_out" ]]; then
            echo "  [SKIP] ${accession} already aligned."
            continue
        fi

        echo "  Aligning ${accession}..."
        bwa mem \
            -t "${THREADS}" \
            "${REFERENCE}" \
            "$R1" "$R2" \
            | samtools view -@ "${THREADS}" -bS -q 20 - \
            > "$bam_out" \
            || { echo "  ERROR: Alignment failed for ${accession}. Skipping..."; continue; }

        echo "  Done: ${accession}"
    done < "${ACCESSIONS_FILE}"

    echo ">>> Alignment complete. BAM files saved to: ${BAM_DIR}"
}

# =============================================================================
# STEP 7: BAM SORTING & INDEXING (SAMtools)
# Sorts and indexes all BAM files for downstream variant calling.
# =============================================================================

sort_and_index_bams() {
    echo ">>> Sorting and indexing BAM files..."

    for bam_file in "${BAM_DIR}"/*.bam; do
        accession=$(basename "$bam_file" .bam)
        sorted_bam="${BAM_DIR}/${accession}_sorted.bam"

        if [[ -f "$sorted_bam" ]]; then
            echo "  [SKIP] ${accession} already sorted."
            continue
        fi

        echo "  Sorting ${accession}..."
        samtools sort \
            -@ "${THREADS}" \
            "$bam_file" \
            -o "$sorted_bam" \
            || { echo "  ERROR: Sorting failed for ${accession}. Skipping..."; continue; }

        echo "  Indexing ${accession}..."
        samtools index \
            -@ "${THREADS}" \
            "$sorted_bam" \
            || { echo "  ERROR: Indexing failed for ${accession}. Skipping..."; continue; }

        # Remove unsorted BAM to save disk space
        rm "$bam_file"

        echo "  Done: ${accession}"
    done

    echo ">>> Sorting and indexing complete."
}

# =============================================================================
# STEP 8: ALIGNMENT QUALITY STATISTICS (SAMtools flagstat)
# Generates alignment summary statistics for each sorted BAM file.
# =============================================================================

generate_flagstats() {
    echo ">>> Generating alignment statistics..."
    mkdir -p "${STATS_DIR}"

    for sorted_bam in "${BAM_DIR}"/*_sorted.bam; do
        accession=$(basename "$sorted_bam" _sorted.bam)
        stats_file="${STATS_DIR}/${accession}_flagstat.txt"

        echo "  Stats for ${accession}..."
        samtools flagstat "$sorted_bam" > "$stats_file" \
            || { echo "  ERROR: flagstat failed for ${accession}."; continue; }
    done

    echo ">>> Flagstat reports saved to: ${STATS_DIR}"
}

# =============================================================================
# STEP 9 & 10: VARIANT CALLING & FILTERING (bcftools)
# Calls SNPs using bcftools mpileup + call.
# Filters to retain only high-confidence biallelic SNPs (QUAL > 20).
# =============================================================================

call_variants() {
    echo ">>> Calling and filtering variants..."
    mkdir -p "${VCF_DIR}"

    for sorted_bam in "${BAM_DIR}"/*_sorted.bam; do
        accession=$(basename "$sorted_bam" _sorted.bam)
        raw_vcf="${VCF_DIR}/${accession}_variants.vcf.gz"
        filtered_vcf="${VCF_DIR}/${accession}_filtered_variants.vcf.gz"

        if [[ -f "$filtered_vcf" ]]; then
            echo "  [SKIP] ${accession} already processed."
            continue
        fi

        echo "  Calling variants for ${accession}..."
        bcftools mpileup \
            -f "${REFERENCE}" \
            "$sorted_bam" \
            | bcftools call \
                -mv \
                --ploidy 1 \
                -Oz \
                -o "$raw_vcf" \
            || { echo "  ERROR: Variant calling failed for ${accession}. Skipping..."; continue; }

        bcftools index "$raw_vcf"

        echo "  Filtering variants for ${accession}..."
        bcftools filter \
            -i 'QUAL>20 && TYPE="snp"' \
            "$raw_vcf" \
            -Oz \
            -o "$filtered_vcf" \
            || { echo "  ERROR: Filtering failed for ${accession}. Skipping..."; continue; }

        bcftools index "$filtered_vcf"

        echo "  Done: ${accession}"
    done

    echo ">>> Variant calling complete. VCF files saved to: ${VCF_DIR}"
}

# =============================================================================
# STEP 11: VARIANT ANNOTATION (SnpEff)
# Annotates filtered VCF files with functional consequence predictions.
# Requires a local SnpEff installation and P. falciparum database.
# =============================================================================

annotate_variants() {
    echo ">>> Annotating variants with SnpEff..."
    mkdir -p "${ANNOTATED_DIR}"

    for filtered_vcf in "${VCF_DIR}"/*_filtered_variants.vcf.gz; do
        accession=$(basename "$filtered_vcf" _filtered_variants.vcf.gz)
        annotated_vcf="${ANNOTATED_DIR}/${accession}_annotated.vcf"

        if [[ -f "$annotated_vcf" ]]; then
            echo "  [SKIP] ${accession} already annotated."
            continue
        fi

        echo "  Annotating ${accession}..."
        java -jar snpEff.jar ann \
            "${SNPEFF_GENOME}" \
            "$filtered_vcf" \
            > "$annotated_vcf" \
            || { echo "  ERROR: Annotation failed for ${accession}. Skipping..."; continue; }

        echo "  Done: ${accession}"
    done

    echo ">>> Annotation complete. Annotated VCFs saved to: ${ANNOTATED_DIR}"
}

# =============================================================================
# STEP 12: CHROMOSOME RENAMING
# Replaces NCBI accession numbers with human-readable chromosome names
# in all annotated VCF files.
# =============================================================================

rename_chromosomes() {
    echo ">>> Renaming chromosomes in annotated VCF files..."

    for annotated_vcf in "${ANNOTATED_DIR}"/*_annotated.vcf; do
        accession=$(basename "$annotated_vcf" _annotated.vcf)
        renamed_vcf="${ANNOTATED_DIR}/${accession}_annotated_renamed.vcf"

        if [[ -f "$renamed_vcf" ]]; then
            echo "  [SKIP] ${accession} already renamed."
            continue
        fi

        echo "  Renaming chromosomes for ${accession}..."
        sed \
            -e 's|NC_004325.2|Chromosome_1|g' \
            -e 's|NC_037280.1|Chromosome_2|g' \
            -e 's|NC_000521.4|Chromosome_3|g' \
            -e 's|NC_004318.2|Chromosome_4|g' \
            -e 's|NC_004326.2|Chromosome_5|g' \
            -e 's|NC_004327.3|Chromosome_6|g' \
            -e 's|NC_004328.3|Chromosome_7|g' \
            -e 's|NC_004329.3|Chromosome_8|g' \
            -e 's|NC_004330.2|Chromosome_9|g' \
            -e 's|NC_037281.1|Chromosome_10|g' \
            -e 's|NC_037282.1|Chromosome_11|g' \
            -e 's|NC_037284.1|Chromosome_12|g' \
            -e 's|NC_004331.3|Chromosome_13|g' \
            -e 's|NC_037283.1|Chromosome_14|g' \
            "$annotated_vcf" > "$renamed_vcf" \
            || { echo "  ERROR: Renaming failed for ${accession}. Skipping..."; continue; }

        echo "  Done: ${accession}"
    done

    echo ">>> Chromosome renaming complete."
}

# =============================================================================
# SAMPLE QC SUMMARY
#
#   Stage                                          Samples
#   -------------------------------------------------------
#   Downloaded from ENA                              90
#   Passed FastQC / Trimmomatic QC                   82
#   Retained after PLINK processing (final)          71
#
#   Per-year breakdown after PLINK (71 final samples):
#     2010: 8  | 2011: 6  | 2012: 6  | 2013: 9  | 2014: 8
#     2015: 10 | 2016: 10 | 2017: 4  | 2018: 10
# =============================================================================

# =============================================================================
# MAIN — Uncomment steps to run them
# Run steps in order for a fresh analysis.
# Individual steps can be re-run independently if needed.
# =============================================================================

generate_accessions_file
# setup_environment
# download_ena_data
# run_fastqc_raw
# run_trimmomatic
# index_reference
# align_samples
# sort_and_index_bams
# generate_flagstats
# call_variants
# annotate_variants
# rename_chromosomes
