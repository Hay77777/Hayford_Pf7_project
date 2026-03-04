# Hayford_Pf7_project

## Longitudinal Genomic Surveillance of Antimalarial Drug Resistance–Associated Single Nucleotide Polymorphisms in *Plasmodium falciparum*, in the Upper East Region of Ghana

**Authors:** Adjoa Agyemang Boakye, Hayford Osei Offei, David Adedia, Enoch Aninagyei

**Affiliation:** Department of Biomedical Sciences, School of Basic and Biomedical Sciences, University of Health and Allied Sciences, Ho, Ghana

**Corresponding Author:** Enoch Aninagyei — eaninagyei@uhas.edu.gh

---

## 📌 Overview

This repository contains the analysis scripts and output figures for a retrospective in silico genomic study of *Plasmodium falciparum* drug resistance in the **Upper East Region of Ghana (2009–2018)**, using data from the [MalariaGEN Pf7 dataset](https://www.malariagen.net/data/pf7) and whole-genome sequencing reads from the [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena).

The study characterises:

- **Temporal allele frequency trends** in resistance-associated SNPs across *pfcrt*, *pfmdr1*, *pfdhfr*, *pfdhps*, and *pfk13*
- **Population diversity and structure** using PCA and neighbor-joining phylogenetic reconstruction
- **Whole-genome sequencing (WGS) analysis** of 90 genomes (2010–2018) including QC, alignment, variant calling, and filtering

**Key findings:**
- The *pfcrt* K76T (CVIET) allele declined steadily and was absent by 2018, indicating near-complete reversion to chloroquine sensitivity
- Antifolate-resistant *pfdhfr/pfdhps* haplotypes persisted throughout the study period
- No validated *pfk13* mutations associated with artemisinin resistance were detected
- Phylogenetic analysis revealed a shift from stable early lineages to a genetically diverse, expanding population post-2012

---

## 🗂️ Repository Structure

```
Hayford_Pf7_project/
├── data/
│   └── README.md                               # Download instructions for input data
├── results/
│   ├── Phylogenetic_Tree_4.pdf                 # Phylogenetic tree (PDF)
│   ├── Phylogenetic_Tree_4_page-0001.jpg       # Phylogenetic tree (image)
│   ├── crt_72.76.CVMNK_allele_frequency.png
│   ├── crt_76.K_allele_frequency.png
│   ├── crt_CVIET_CVMNK_trend_plot_fixed.png    # Manuscript figure
│   ├── dhfr_108.S_allele_frequency.png
│   ├── dhfr_164.I_allele_frequency.png
│   ├── dhfr_51.N_allele_frequency.png
│   ├── dhfr_59.C_allele_frequency.png
│   ├── dhps_437.G_allele_frequency.png
│   ├── dhps_540.K_allele_frequency.png
│   ├── dhps_581.A_allele_frequency.png
│   ├── dhps_613.A_allele_frequency.png
│   ├── double_triple_mutant_plot_clean.png      # Manuscript figure
│   ├── k13_mutation_trend_plot_fixed.png        # Manuscript figure
│   ├── kelch13_349.726_ns_changes_allele_frequency.png
│   ├── mdr1_breakpoint_allele_frequency.png
│   ├── mdr1_dup_call_allele_frequency.png
│   ├── partially_resistant_haplotype_plot2.png  # Manuscript figure
│   └── pm2_breakpoint_allele_frequency.png
├── scripts/
│   ├── plasmodium_variant_pipeline.sh          # Steps 1–12: download → annotation
│   ├── 04_filter_and_summarize_snps.sh         # SNP filtering by resistance gene
│   ├── 01_plink_pca_analysis.R                 # PLINK QC + PCA + clustering
│   ├── 02_phylogenetic_tree.R                  # Neighbor-joining phylogenetic tree
│   └── 03_snp_trend_analysis.R                 # SNP allele frequency trend analysis
└── README.md
```

---

## 📥 Input Data

### MalariaGEN Pf7 (SNP trend analysis — 438 samples)

| File | Description |
|------|-------------|
| `Pf7_samples.txt` | Sample metadata including country, region, and year of collection |
| `Pf7_drug_resistance_marker_genotypes.txt` | Pre-curated genotype calls at key drug resistance loci |

> **Download:** https://www.malariagen.net/data/pf7
> See `data/README.md` for full instructions.

### European Nucleotide Archive (WGS analysis — 90 samples)

Raw FASTQ reads for the WGS subset were downloaded from ENA. Sample accession IDs correspond to MalariaGEN Pf7 Ghana / Upper East Region isolates (2010–2018).

> **ENA Portal:** https://www.ebi.ac.uk/ena
> See `data/README.md` for the full list of 90 ENA accession numbers.

Samples were filtered to retain only those from the **Kassena Nankana West District, Upper East Region, Ghana**.

---

## 🔬 Analysis Pipeline

### 1. WGS Pipeline (`plasmodium_variant_pipeline.sh`)
End-to-end shell pipeline covering:

| Step | Tool | Details |
|------|------|---------|
| Download | wget (ENA FTP) | Paired-end FASTQ for all 90 ERR samples |
| QC | FastQC | Per-base quality, GC content, adapter contamination |
| Trimming | Trimmomatic | Sliding-window 4:20, min Q20 |
| Alignment | BWA-MEM v0.7.17 | Mapped to *Pf* 3D7 reference (GCF_000002765.6); MAPQ ≥ Q20 |
| BAM processing | SAMtools v1.17 | Sort, index, flagstat QC |
| Variant calling | bcftools v1.17 | mpileup + call; QUAL > 20; biallelic SNPs only |
| Annotation | SnpEff | Functional consequence prediction |
| Chr renaming | sed | NCBI accessions → chromosome names |

**Sample QC summary:**

| Stage | Samples |
|-------|---------|
| Downloaded from ENA | 90 |
| Passed FastQC / Trimmomatic | 82 |
| Retained after PLINK processing | 71 |

### 2. SNP Filtering (`04_filter_and_summarize_snps.sh`)
Filters SNPs for 7 resistance genes (*pfcrt*, *pfmdr1*, *pfdhps*, *pfdhfr*, *pfpm2*, *pfk13*) and generates summary statistics per gene, per year, and amino acid change frequencies.

### 3. PLINK QC & PCA (`01_plink_pca_analysis.R`)
- Converts SNP data to PLINK binary format
- Imputes missing genotypes (`--fill-missing-a2`)
- Runs PCA (10 components) on 71 retained samples and 69 SNPs
- Hierarchical clustering with Ward's method
- ANOVA and cophenetic correlation statistics

### 4. Phylogenetic Tree (`02_phylogenetic_tree.R`)
- Builds binary SNP presence/absence matrix
- Computes Jaccard genetic distances
- Constructs neighbor-joining tree using `ape`
- Tips coloured by collection year

### 5. SNP Trend Analysis (`03_snp_trend_analysis.R`)
- Merges MalariaGEN metadata with per-year filtered samples
- Generates allele frequency trend plots for all SNP markers
- Manuscript figures: *pfcrt* CVIET/CVMNK, *pfdhfr* double/triple mutant haplotypes
- Pearson's chi-square trend tests (p-values)

---

## 📊 Manuscript Figures

| Figure | File |
|--------|------|
| *pfcrt* CVIET/CVMNK haplotype trends | `results/crt_CVIET_CVMNK_trend_plot_fixed.png` |
| *pfdhfr* double & triple mutant haplotypes | `results/double_triple_mutant_plot_clean.png` |
| *pfk13* mutation trends | `results/k13_mutation_trend_plot_fixed.png` |
| Partially resistant haplotype frequencies | `results/partially_resistant_haplotype_plot2.png` |
| Phylogenetic tree | `results/Phylogenetic_Tree_4.pdf` |

---

## 🛠️ Dependencies

### R packages
| Package | Purpose |
|---------|---------|
| R v4.1.2 | Core analysis environment |
| ggplot2 | Visualisation |
| ggrepel | Outlier label repelling in PCA |
| dplyr / tidyr | Data manipulation |
| ape | Phylogenetic tree construction |
| ggdendro / cluster | Hierarchical clustering |

### Bioinformatics tools
| Tool | Version |
|------|---------|
| FastQC | — |
| Trimmomatic | — |
| BWA-MEM | v0.7.17-r1188 |
| SAMtools | v1.17 |
| bcftools | v1.17 |
| PLINK | v1.9 |
| SnpEff | — |

> All analyses were performed in a **WSL2 Ubuntu** environment on an Intel Core i7 system with 16 GB RAM.

---

## ▶️ How to Reproduce

```bash
# 1. Clone the repository
git clone https://github.com/Hay77777/Hayford_Pf7_project.git
cd Hayford_Pf7_project

# 2. Download MalariaGEN Pf7 input files and place in data/
#    See data/README.md for full instructions

# 3. Run the WGS pipeline (Steps 1-12)
bash scripts/plasmodium_variant_pipeline.sh

# 4. Filter SNPs by resistance gene
bash scripts/04_filter_and_summarize_snps.sh

# 5. Run R analyses (in order)
Rscript scripts/01_plink_pca_analysis.R
Rscript scripts/02_phylogenetic_tree.R
Rscript scripts/03_snp_trend_analysis.R
```

---

## 📄 Citation

If you use this code or results, please cite:

> Agyemang Boakye A, Osei Offei H, Adecia D, Aninagyei E. (*under review*). Longitudinal Genomic Surveillance of Antimalarial Drug Resistance–Associated Single Nucleotide Polymorphisms in *Plasmodium falciparum*, in the Upper East Region of Ghana.

> MalariaGEN *et al.* (2023). An open dataset of *Plasmodium falciparum* genome variation in 7,000 worldwide samples. *Wellcome Open Research*. https://doi.org/10.12688/wellcomeopenres.18681.1

---

## 🔗 Data Sources

| Resource | Link |
|----------|------|
| MalariaGEN Pf7 | https://www.malariagen.net/data/pf7 |
| Pf7 DOI | https://doi.org/10.12688/wellcomeopenres.18681.1 |
| ENA | https://www.ebi.ac.uk/ena |
| *Pf* 3D7 Reference Genome | GCF_000002765.6 |

---


