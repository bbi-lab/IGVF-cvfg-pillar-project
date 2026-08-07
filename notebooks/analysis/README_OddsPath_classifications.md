# OddsPath Classifications Pipeline

This Jupyter notebook performs **ACMG/AMP-style variant classification** using **OddsPath functional evidence** combined with **genome-wide computational predictor calibrations** (REVEL, AlphaMissense, MutPred2).

## Overview

This pipeline is a companion to the `Variant_Classification_analysis.ipynb` pipeline, but uses OddsPath-derived functional points (`OP_points`) instead of ExCALIBR calibrations. This approach is specifically designed for datasets where OddsPath likelihood ratios have been calculated.

The pipeline:
1. Loads the integrated variant effect analysis dataset
2. Identifies and flags conflicting functional evidence using OddsPath points
3. Deduplicates variants by selecting those with maximum absolute OddsPath points
4. Combines OddsPath functional evidence with genome-wide predictor scores
5. Generates classifications for different variant categories

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
| `integrated_variant_effect_dataset_analysis.csv.gz` | `outputs/` | Output from the Variant Classification Analysis pipeline |

## Methods

### Key Differences from Variant_Classification_analysis Pipeline

| Aspect | Variant_Classification_analysis | OddsPath_classifications |
|--------|--------------------------------|--------------------------|
| Functional evidence | ExCALIBR calibrations (`Fxn_points`) | OddsPath LRs (`OP_points`) |
| Predictor calibrations | Gene-specific (with genome-wide fallback) | Genome-wide only |
| Conflict detection | Based on `Fxn_points` | Based on `OP_points` |
| Applicable datasets | All datasets with ExCALIBR calibrations | Datasets with OddsPath calibrations |

### OddsPath Points

OddsPath likelihood ratios are converted to evidence points:

| OddsPath (LR) | Direction | Points | Evidence Code |
|---------------|-----------|--------|---------------|
| > 350 | Abnormal | +8 | PS3_very_strong |
| > 18.7 | Abnormal | +4 | PS3_strong |
| > 4.3 | Abnormal | +2 | PS3_moderate |
| > 2.1 | Abnormal | +1 | PS3_supporting |
| < 0.48 | Normal | −1 | BS3_supporting |
| < 0.23 | Normal | −2 | BS3_moderate |
| < 0.053 | Normal | −4 | BS3_strong |

### Classification Thresholds

Total points (OddsPath + Predictor) are mapped to classifications:

| Total Points | Classification |
|--------------|----------------|
| ≥ 10 | Pathogenic |
| 6 – 9 | Likely Pathogenic |
| 0 – 5 | Uncertain Significance |
| −6 – −1 | Likely Benign |
| ≤ −7 | Benign |

### Variant Deduplication

When the same variant appears in multiple assays, the pipeline selects a representative using the **maximum absolute OddsPath points** strategy:

1. **Group variants** by genomic coordinates (nucleotide-level) or amino acid position (protein-level)
2. **Detect conflicting evidence** — Flag variants where different assays give opposite-direction OddsPath points
3. **Select representative** — Keep the variant with the highest absolute `OP_points` value
4. **Apply assay priority** — For tied OddsPath points in amino acid assays, use predefined assay priority list

### Excluded Variants

Before deduplication, the following are removed:
- SFPQ gene (insufficient controls)
- Splice variants
- Start-lost variants in amino acid assays
- ClinVar conflict flags
- Predictor training variants (for REVEL and MutPred2 analyses)

## Usage

### 1. Set Project Root

```bash
export PROJECT_ROOT=/path/to/your/project
jupyter notebook OddsPath_classifications.ipynb
```

### 2. Prerequisites

Ensure these pipelines have been run first:
1. `Integrated_variant_effect_dataset_pipeline.ipynb`
2. `OddsPath_calculations.ipynb`
3. `Variant_Classification_analysis.ipynb`

### 3. Execute the Pipeline

Run cells sequentially. The notebook will:
1. Load the integrated variant effect analysis dataset
2. Flag splice variants and start-lost variants
3. Detect conflicting OddsPath evidence between assays
4. Select variants with maximum absolute OddsPath points
5. Subset variants by category (controls, VUS, gnomAD, ClinGen, unobserved)
6. Export per-category classification files

## Output

### Output Files

```
outputs/OddsPath_calibrations/
├── controls_REVEL_OP.csv
├── controls_MP2_OP.csv
├── controls_AM_OP.csv
├── VUS_REVEL_OP.csv
├── VUS_MP2_OP.csv
├── VUS_AM_OP.csv
├── ClinGen_Repo_REVEL_OP.csv
├── ClinGen_repo_MP2_OP.csv
├── ClinGen_repo_AM_OP.csv
├── gnomAD_REVEL_OP.csv
├── gnomAD_AM_OP.csv
├── gnomAD_MP2_OP.csv
├── Unobserved_REVEL_OP.csv
├── Unobserved_AM_OP.csv
└── Unobserved_mut_OP.csv

outputs/Supplementary_Data_6.xlsx.gz  # All category files combined into Excel
```

### Key Output Columns
 See column descriptions for Supp Dat 6

## Variant Categories

| Category | Description | Filter Criteria |
|----------|-------------|-----------------|
| **controls** | ClinVar P/LP and B/LB variants | `clnsig_group_18_25` ∈ {Pathogenic, Likely pathogenic, Benign, Likely benign} with 1+ stars |
| **VUS** | Variants of Uncertain Significance | `clinvar_18_25` = 'Uncertain significance' |
| **ClinGen_Repo** | ClinGen Evidence Repository curated variants | `Updated_Classification_ClinGen_repo` is not null and not VUS |
| **gnomAD** | Population variants from gnomAD | `gnomad_MAF` is not null |
| **Unobserved** | Variants not seen in ClinVar or gnomAD | `clinvar_sig_2025` is null AND `gnomad_MAF` is null (SNVs only) |

## Notes

- This pipeline uses **genome-wide** predictor calibrations only (no gene-specific thresholds)
- OddsPath points are derived from likelihood ratios calculated in the `OddsPath_calculations.ipynb` pipeline
- The `StandardizedClass` column determines whether `OddsNormal` or `OddsAbnormal` is used for point assignment
- Variants with missing OddsPath values receive 0 functional points
- REVEL and MutPred2 training variants are excluded from their respective analyses
- ClinVar significance uses 2018 controls for BRCA1, PTEN, MSH2, TP53; 2025 controls for all other genes
