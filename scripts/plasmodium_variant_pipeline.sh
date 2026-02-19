#!/bin/bash
# =============================================================================
# Plasmodium falciparum Variant Calling Pipeline
# Project: PRJNA705555
# Reference Genome: GCF_000002765.6 (P. falciparum 3D7)
# =============================================================================
# PIPELINE OVERVIEW:
#   1. Environment Setup
#   2. Data Download from SRA
#   3. Quality Control (FastQC)
#   4. Alignment to Reference Genome (BWA-MEM)
#   5. SAM to BAM Conversion & Indexing (SAMtools)
#   [6. SNP Annotation - SnpEff (to be added)]
# =============================================================================

set -euo pipefail  # Exit on error, treat unset variables as errors

# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_ID="PRJNA705555"
REFERENCE_GENOME="GCF_000002765.6_GCA_000002765_genomic.fna"
FASTQ_DIR="Fastq_files"
SAM_OUTPUT_DIR="SAM_Output"
FASTQC_OUTPUT_DIR="fastqc_output"

SAMPLE_IDS=(
    "SRR13822575"
    "SRR13822576"
    "SRR13822577"
    "SRR13822578"
    "SRR13822579"
    "SRR13822580"
    "SRR13822581"
    "SRR13822582"
    "SRR13822583"
)

# =============================================================================
# STEP 1: ENVIRONMENT SETUP
# Install required tools (run once)
# =============================================================================

setup_environment() {
    echo ">>> Setting up environment..."
    sudo apt update
    sudo apt install -y build-essential libncbi-vdb-dev sra-toolkit fastqc

    echo ">>> Verifying SRA Toolkit installation..."
    fastq-dump --version
}

# =============================================================================
# STEP 2: DATA DOWNLOAD FROM SRA
# =============================================================================

download_sra_data() {
    echo ">>> Fetching SRA run info for project: ${PROJECT_ID}..."
    esearch -db sra -query "${PROJECT_ID}" | efetch -format runinfo > runinfo.csv

    echo ">>> Extracting SRA run IDs..."
    cut -d ',' -f 1 runinfo.csv | tail -n +2 > sra_ids.txt

    echo ">>> Downloading .sra files with prefetch..."
    while read -r sra_id; do
        echo "  Downloading ${sra_id}..."
        prefetch "${sra_id}"
    done < sra_ids.txt

    echo ">>> Converting .sra files to FASTQ format..."
    while read -r sra_id; do
        echo "  Converting ${sra_id}..."
        fasterq-dump "${sra_id}" --outdir "${FASTQ_DIR}"
    done < sra_ids.txt

    # Verify downloaded files
    echo ">>> Listing downloaded .sra files:"
    find . -name "*.sra" -exec ls -lh {} \;
}

# =============================================================================
# STEP 3: QUALITY CONTROL
# Reference: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000002765.6/
# Note: SRR13822576 flagged for failed per-sequence GC content and 
#       sequence duplication levels. Retained for downstream comparison.
# =============================================================================

run_fastqc() {
    echo ">>> Running FastQC on all FASTQ files..."
    mkdir -p "${FASTQC_OUTPUT_DIR}"
    fastqc "${FASTQ_DIR}"/*.fastq --outdir "${FASTQC_OUTPUT_DIR}"
    echo ">>> FastQC reports saved to: ${FASTQC_OUTPUT_DIR}"
}

# =============================================================================
# STEP 4 & 5: ALIGNMENT & BAM CONVERSION
# BWA-MEM alignment followed by SAMtools conversion and indexing
# =============================================================================

index_reference() {
    echo ">>> Indexing reference genome: ${REFERENCE_GENOME}..."
    bwa index "${REFERENCE_GENOME}"
}

align_samples() {
    echo ">>> Starting alignment for all samples..."
    mkdir -p "${SAM_OUTPUT_DIR}"

    for sample_id in "${SAMPLE_IDS[@]}"; do
        echo "  Processing ${sample_id}..."

        forward_read="${FASTQ_DIR}/${sample_id}_1.fastq"
        reverse_read="${FASTQ_DIR}/${sample_id}_2.fastq"
        bam_file="${SAM_OUTPUT_DIR}/${sample_id}.bam"
        tmp_sam="${SAM_OUTPUT_DIR}/${sample_id}.sam"

        # Validate input files
        if [[ ! -f "${forward_read}" || ! -f "${reverse_read}" ]]; then
            echo "  WARNING: Missing FASTQ files for ${sample_id}. Skipping..."
            continue
        fi

        # Align with BWA-MEM
        bwa mem "${REFERENCE_GENOME}" "${forward_read}" "${reverse_read}" > "${tmp_sam}" \
            || { echo "  ERROR: Alignment failed for ${sample_id}"; continue; }

        # Convert SAM to BAM
        samtools view -bS "${tmp_sam}" > "${bam_file}" \
            || { echo "  ERROR: SAM to BAM conversion failed for ${sample_id}"; continue; }

        # Index BAM file
        samtools index "${bam_file}" \
            || { echo "  ERROR: BAM indexing failed for ${sample_id}"; continue; }

        # Remove intermediate SAM file to save disk space
        rm "${tmp_sam}"

        echo "  Done: ${sample_id}"
    done

    echo ">>> All samples processed. BAM files saved to: ${SAM_OUTPUT_DIR}"
}

# =============================================================================
# MAIN — Comment/uncomment steps as needed
# =============================================================================

# setup_environment
# download_sra_data
# run_fastqc
# index_reference
# align_samples
