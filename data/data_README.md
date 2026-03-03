# Data

This directory contains **no raw data files**. Input data must be downloaded directly from MalariaGEN and the European Nucleotide Archive (ENA).

---

## 📥 Required Input Files

This analysis uses data from two public sources:

| Source | Files | Used for |
|--------|-------|---------|
| MalariaGEN Pf7 | `Pf7_samples.txt`, `Pf7_drug_resistance_marker_genotypes.txt` | SNP trend analysis (438 samples) |
| ENA | FASTQ files (ERR accessions below) | WGS analysis (90 samples) |

---

## 1. MalariaGEN Pf7 Files

### How to Download
1. Visit: **https://www.malariagen.net/data/pf7**
2. Navigate to the **"Data downloads"** section
3. Download both:
   - `Pf7_samples.txt`
   - `Pf7_drug_resistance_marker_genotypes.txt`
4. Place both files in this `data/` directory

### Citation
> MalariaGEN et al. (2023). An open dataset of *Plasmodium falciparum* genome variation
> in 7,000 worldwide samples. *Wellcome Open Research*.
> https://doi.org/10.12688/wellcomeopenres.18681.1

---

## 2. WGS Samples — ENA Accession Numbers

90 whole-genome sequencing samples from the **Kassena Nankana West District, Upper East Region, Ghana (2010–2018)** were downloaded from the European Nucleotide Archive (ENA).

**ENA Portal:** https://www.ebi.ac.uk/ena

Search each accession number on the ENA portal and download the paired-end FASTQ files.

| Year | ENA Accession |
|------|--------------|
| 2010 | ERR019549 |
| 2010 | ERR039235 |
| 2010 | ERR036147 |
| 2010 | ERR035362 |
| 2010 | ERR045604 |
| 2010 | ERR205933 |
| 2010 | ERR211562 |
| 2010 | ERR045620 |
| 2010 | ERR045607 |
| 2010 | ERR045609 |
| 2011 | ERR039979 |
| 2011 | ERR042698 |
| 2011 | ERR042701 |
| 2011 | ERR042719 |
| 2011 | ERR042722 |
| 2011 | ERR042727 |
| 2011 | ERR063567 |
| 2011 | ERR114406 |
| 2011 | ERR114409 |
| 2011 | ERR211462 |
| 2012 | ERR376187 |
| 2012 | ERR376188 |
| 2012 | ERR376201 |
| 2012 | ERR376202 |
| 2012 | ERR376203 |
| 2012 | ERR376212 |
| 2012 | ERR376213 |
| 2012 | ERR376214 |
| 2012 | ERR376216 |
| 2012 | ERR403196 |
| 2013 | ERR450116 |
| 2013 | ERR450108 |
| 2013 | ERR450099 |
| 2013 | ERR450045 |
| 2013 | ERR586212 |
| 2013 | ERR586179 |
| 2013 | ERR450104 |
| 2013 | ERR450112 |
| 2013 | ERR450044 |
| 2013 | ERR450041 |
| 2014 | ERR3504126 |
| 2014 | ERR3523699 |
| 2014 | ERR3546199 |
| 2014 | ERR3546151 |
| 2014 | ERR3546188 |
| 2014 | ERR3594524 |
| 2014 | ERR3594529 |
| 2014 | ERR3523718 |
| 2014 | ERR3523721 |
| 2014 | ERR3594576 |
| 2015 | ERR3486201 |
| 2015 | ERR3523571 |
| 2015 | ERR3523563 |
| 2015 | ERR3546220 |
| 2015 | ERR3523635 |
| 2015 | ERR3523614 |
| 2015 | ERR1214163 |
| 2015 | ERR1214188 |
| 2015 | ERR1214209 |
| 2015 | ERR1214148 |
| 2016 | ERR2532525 |
| 2016 | ERR2542018 |
| 2016 | ERR2541936 |
| 2016 | ERR2532503 |
| 2016 | ERR2532517 |
| 2016 | ERR2532524 |
| 2016 | ERR2532570 |
| 2016 | ERR2541897 |
| 2016 | ERR2541942 |
| 2016 | ERR2532534 |
| 2017 | ERR3486300 |
| 2017 | ERR2891366 |
| 2017 | ERR2890914 |
| 2017 | ERR3608767 |
| 2017 | ERR3546081 |
| 2017 | ERR2375226 |
| 2017 | ERR2506008 |
| 2017 | ERR3504084 |
| 2017 | ERR3523658 |
| 2017 | ERR3504120 |
| 2018 | ERR3608800 |
| 2018 | ERR3608865 |
| 2018 | ERR3584487 |
| 2018 | ERR3594738 |
| 2018 | ERR3608785 |
| 2018 | ERR3594774 |
| 2018 | ERR3594757 |
| 2018 | ERR3594689 |
| 2018 | ERR3594679 |
| 2018 | ERR3594661 |

---

## ✅ After Downloading

Your `data/` folder should look like this before running the pipeline:

```
data/
├── README.md
├── Pf7_samples.txt
├── Pf7_drug_resistance_marker_genotypes.txt
└── fastq/
    ├── ERR019549_1.fastq.gz
    ├── ERR019549_2.fastq.gz
    └── ... (paired files for all 90 samples)
```

Then run the pipeline from the project root:

```bash
bash scripts/plasmodium_variant_pipeline.sh
```
