# Variant Classification Analysis Pipeline

This Jupyter notebook performs **ACMG/AMP-style variant classification** by combining functional assay evidence (PS3/BS3) with computational predictor evidence (PP3/BP4) using a Bayesian point-based framework.

## Overview

The pipeline integrates multiple evidence sources to classify variants:

1. **Functional evidence** from MAVE assays (via ExCALIBR calibrations or OddsPath (for F9 and TP53))
2. **Computational predictor evidence** from REVEL, AlphaMissense, and MutPred2
3. Both **genome-wide** and **gene-specific** calibrations for predictors

Evidence is converted to points and summed to produce final classifications following the ACMG/AMP framework.

The pipeline also performs **deduplication** by selecting variants with the highest absolute point values and generates separate output files for different variant categories: ClinVar controls, ClinGen-curated variants, VUS, gnomAD variants, and unobserved variants.

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

## Input Data

### Required Files

| File | Location | Description |
|------|----------|-------------|
| `integrated_variant_effect_dataset.tsv.gz` | `outputs/` | Output from the Integrated Variant Effect Dataset pipeline |
| `Supplementary_Data_3.xlsx` | `data/Supplemental_Data/` | Dataset curation sheet |
| `Supplementary_Data_4.xlsx` | `data/Supplemental_Data/` | Calibration data (multiple sheets) |

### Supplementary Data 4 Sheets

| Sheet Name | Description |
|------------|-------------|
| `ExCALIBR_calibrations` | Functional assay score interval calibrations |
| `REVEL_gene_specific_calibration` | Gene-specific REVEL thresholds |
| `AM_gene_specific_calibrations` | Gene-specific AlphaMissense thresholds |
| `MP2_gene_specific_calibrations` | Gene-specific MutPred2 thresholds |
| `OddsPath_calibrations` | OddsPath likelihood ratios per dataset |

## Methods

### Evidence Point System

Evidence codes are converted to points following Tavtigian et al. 2020 and for predictors thresholds from Bergquist et al. are used

### Gene-Specific vs Genome-Wide Calibrations

The pipeline preferentially uses **gene-specific** predictor thresholds when available, falling back to genome-wide thresholds otherwise.

### ClinVar Control Selection

- **2018 ClinVar controls**: Used for BRCA1, PTEN, MSH2 (legacy calibrations)
- **2025 ClinVar controls**: Used for all other genes


## Variant Deduplication

When the same genomic variant appears in multiple assays, the pipeline selects a single representative using the **maximum absolute functional points** strategy:

### Deduplication Logic

1. **Group variants** by genomic coordinates:
   - Nucleotide-level assays: `Gene`, `Chrom`, `hg38_start`, `ref_allele`, `alt_allele`
   - Protein-level assays: `Gene`, `aa_pos`, `aa_ref`, `aa_alt`, `RefSeq_transcript_ID`

2. **Identify conflicting functional data** — Variants where different assays give opposite-direction evidence (one positive, one negative) are flagged as `conflicting_fxn_data`

3. **Select representative variant**:
   - For non-conflicting variants: Keep the row with the **highest absolute `Fxn_points`** value
   - Mark selected rows as `First_max_fxn_pts` (nucleotide) or `max_fxn_pts` (protein-level)
   - Similarly select maximum predictor points per variant (`max_pred_pts`)

### Conflict Detection

The pipeline flags cases where functional and computational evidence conflict (point in opposite directions). Retains the variant, but flags it. The `Conflicting_*` columns indicate the combined point total when evidence conflicts.

### Excluded Variants

Before deduplication, the following variants are removed:
- Flagged variants (`Flag = '*'`)
- Splice variants
- Conflicting functional data
- Start-lost variants in amino acid assays
- Predictor training variants (for REVEL and MutPred2 analyses)

## Variant Category Outputs

The pipeline generates separate deduplicated output files for each variant category and predictor combination:

### Output Categories

| Category | Description | Filter Criteria |
|----------|-------------|-----------------|
| **controls** | ClinVar P/LP and B/LB variants | `clnsig_group_18_25` ∈ {Pathogenic, Likely pathogenic, Benign, Likely benign} with 1+ stars |
| **ClinGen_Repo** | ClinGen Evidence Repository curated variants | `Updated_Classification_ClinGen_repo` is not null and not VUS |
| **VUS** | Variants of Uncertain Significance | `clinvar_sig_2025` = 'Uncertain significance' |
| **gnomAD** | Population variants from gnomAD | `gnomad_MAF` is not null |
| **Unobserved** | Variants not seen in ClinVar or gnomAD | `clinvar_sig_2025` is null AND `gnomad_MAF` is null (SNVs only) |

### Per-Predictor Files

Each category is exported separately for each predictor to exclude training variants:

```
outputs/GeneSpecific_calibrations/
├── controls_REVEL_GeneSpecific.csv
├── controls_AM_GeneSpecific.csv
├── controls_MP2_GeneSpecific.csv
├── ClinGen_Repo_REVEL_GeneSpecific.csv
├── ClinGen_Repo_AM_GeneSpecific.csv
├── ClinGen_Repo_MP2_GeneSpecific.csv
├── VUS_REVEL.csv
├── VUS_AM.csv
├── VUS_MP2.csv
├── gnomAD_REVEL.csv
├── gnomAD_AM.csv
├── gnomAD_MP2.csv
├── Unobserved_REVEL.csv
├── Unobserved_AM.csv
└── Unobserved_MP2.csv
```


## Usage

### 1. Set Project Root

```bash
export PROJECT_ROOT=/path/to/your/project
jupyter notebook Variant_Classification_analysis.ipynb
```

### 2. Prerequisites

Ensure these pipelines have been run first:
1. `Integrated_variant_effect_dataset_pipeline.ipynb`
2. `OddsPath_calculations.ipynb`

### 3. Execute the Pipeline

Run cells sequentially. The notebook will:
1. Load integrated variant data and calibration tables
2. Apply ExCALIBR functional score calibrations
3. Merge OddsPath likelihood ratios
4. Standardize functional classifications
5. Apply genome-wide predictor thresholds
6. Apply gene-specific predictor thresholds
7. Calculate total points and assign classifications
8. Detect conflicting evidence
9. Export per-gene classification files

## Output

### Output Files

```
outputs/integrated_variant_effect_dataset_analysis.csv.gz  # Full analysis dataset

outputs/GeneSpecific_calibrations/
├── controls_REVEL_GeneSpecific.csv
├── controls_AM_GeneSpecific.csv
├── controls_MP2_GeneSpecific.csv
├── ClinGen_Repo_REVEL_GeneSpecific.csv
├── ClinGen_Repo_AM_GeneSpecific.csv
├── ClinGen_Repo_MP2_GeneSpecific.csv
├── VUS_REVEL.csv
├── VUS_AM.csv
├── VUS_MP2.csv
├── gnomAD_REVEL.csv
├── gnomAD_AM.csv
├── gnomAD_MP2.csv
├── Unobserved_REVEL.csv
├── Unobserved_AM.csv
└── Unobserved_MP2.csv

outputs/Supplementary_Data_5.xlsx.gz  # All category files combined into Excel
```

### Key Output Columns

See column descriptions for Supp Data 5

### Classification Column Values

The `Class_*` columns contain one of:
- **Pathogenic** — Total points ≥ 10
- **Likely Pathogenic** — Total points 6–9
- **Uncertain** — Total points 0–5
- **Likely Benign** — Total points −6 to −1
- **Benign** — Total points ≤ −7


## Notes

- Variants with missing functional scores receive 0 functional points
- Gene-specific predictor thresholds take precedence over genome-wide when available
- REVEL and MutPred2 training variants are excluded from their respective analyses to prevent circularity
- Unobserved variants are restricted to SNVs only
- TP53 and F9 datasets are limited to meta-analysis versions (individual assay datasets excluded)
- SFPQ gene is excluded due to insufficient controls

