# Integrated Variant Effect Dataset Pipeline

This Jupyter notebook generates **Supplementary Data 1** — an integrated dataset of variant effect measurements from experimental, annotated with genomic coordinates, functional consequence predictions, and computational pathogenicity scores.

## Overview

The pipeline harmonizes variant-level data from 83 published MAVE datasets, standardizes variant nomenclature (HGVS), maps variants to genomic coordinates (GRCh38), and integrates annotations from:

- **Ensembl VEP** (Variant Effect Predictor)
- **VariantValidator** 
- **ClinVar** (clinical significance)
- **gnomAD** (population allele frequencies)
- **Computational predictors**: REVEL, AlphaMissense, MutPred2, SpliceAI

## Requirements

### Software
- Python 3.12+
- Jupyter Notebook/Lab

### Python Dependencies
```
pandas
numpy
openpyxl
```

### External Tools
- **Ensembl VEP** (v113 or compatible) — for variant annotation
- **VariantValidator** — for HGVS validation and coordinate resolution

## Data Availability

### Provided Data (GitHub/Zenodo)

The following input files are hosted in the repository or on Zenodo (https://zenodo.org/records/18637474):

| File/Directory | Description |
|----------------|-------------|
| `data/Variant_score_files/` | Original supplementary tables from published MAVE studies (Excel files), any manual changes are documented in Supp_data_3 in the additional notes column |
| `data/VariantValidator_output/` | Pre-computed VariantValidator results (`variant_validator_output_1.csv.gz`, `variant_validator_output_2.csv.gz`) |
| `data/ClinVar_files/` | ClinVar variant summary files (January 2025 and December 2018|
| `data/gnomAD_files/` | gnomAD allele frequency annotations downloaded form web, v4, only for the 40 genes annotated|
| `data/Simplified_consequence_file/` | Extended Ensembl consequence mapping to simplified consequence (`extended_ensembl_consequence.csv.gz`) |
| `data/Liftover_files/` | h38 <-> hg19 liftover files |
| `data/VEP_outputs/` | Pre-computed VEP annotation files |
| `data/ClinGen_evidence_repo_files/` | ClinGen evidence repository |

### Not Provided (Obtain Separately, too large to host)

Due to file size constraints, the following predictor score files are **not included** and must be obtained from their original sources:

| File | Source | Description |
|------|--------|-------------|
| `REVEL_scores.csv.gz` | REVEL pathogenicity scores |
| `AlphaMissense_hg38.tsv.gz` | AlphaMissense predictions |
| MutPred2 scores | MutPred2 pathogenicity predictions |
| SpliceAI scores | Splice site predictions |

Place these files in `data/Predictor_score_files/` before running the notebook.

## Directory Structure

```
project_root/
├── data/
│   ├── Variant_score_files/                    # MAVE study supplementary files
│   ├── VariantValidator_output/        # VariantValidator results
│   ├── ClinVar_files/                  # ClinVar annotations
│   ├── gnomAD_files/                   # gnomAD allele frequencies
│   ├── VEP_output/                    # VEP annotation outputs
│   ├── Predictor_score_files/          # Pathogenicity predictors (user-supplied)
│   ├── Simplified_consequence_file/    # Consequence term mapping
│   └── ClinGen_evidence_repo_files/    # ClinGen evidence repository files
|   └── Liftover_files/    		# Liftover files
├── tmp/                                # Intermediate files (auto-created)
├── outputs/                            # Final output files (auto-created)
└── Integrated_variant_effect_dataset_pipeline.ipynb
```

## Usage

### 1. Set Project Root

Set the `PROJECT_ROOT` environment variable or run the notebook from the project directory:

```bash
export PROJECT_ROOT=/path/to/your/project
jupyter notebook Integrated_variant_effect_dataset_pipeline.ipynb
```

### 2. Run the Pipeline

Execute cells sequentially. The pipeline will:

1. **Load and harmonize** variant data from 83 MAVE datasets
2. **Standardize nomenclature** (HGVS.p, HGVS.c, amino acid codes)
3. **Map to GRCh38 coordinates** using VEP and VariantValidator
4. **Merge annotations** from ClinVar, gnomAD, and predictors
5. **Export** the integrated dataset

### 3. Output

The final output is written to:

```
outputs/integrated_variant_effect_dataset.tsv.gz
```

## Output Schema

Refer to column descriptions in Supplementary Data 1

## Included Datasets

The pipeline integrates data from 83 MAVE datasets covering genes including:

ASPA, BAP1, BARD1, BRCA1, BRCA2, CALM1/2/3, CARD11, CBS, CHEK2, CRX, CTCF, DDX3X, F9, FKRP, G6PD, GCK, HMBS, JAG1, KCNE1, KCNH2, KCNQ4, LARGE1, MSH2, NDUFAF6, OTC, PALB2, PAX6, PTEN, RAD51C, RAD51D, RHO, SCN5A, SFPQ, SGCB, TARDBP, TP53, TPK1, TSC2, VHL, XRCC2

## Notes

- Variants that cannot be unambiguously mapped are flagged in the `Flag` column, please filter out when using the data
- The `nucleotide_or_aa` column indicates whether the original data was reported at nucleotide or protein level 

