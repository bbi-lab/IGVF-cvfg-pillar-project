
---

##  A scalable approach to resolving variants of uncertain significance

This repository contains all the scripts necessary to recreate the analyses presented in `A scalable approach to resolving variants of uncertain significance`

### 1. Analysis Scripts (`Analysis/`)

#### ▸ Integrated Variant Effect Dataset Pipeline

Scripts used to:

- Harmonize variant identifiers across published experimental datasets  
- Map protein-level variants to nucleotide-level representations    
- Annotate variants with ClinVar, variant effect predictor (REVEL, MutPred2, AlphaMissense) information, and spliceAI scores
-Create an integrated variant effect dataset for XX variants and XX measurements, presented in an expanded format of XX (all possible nucleotide variants mapped to protein variants)

#### ▸ Classification Analysis

Scripts implementing:

- OddsPath calculations for datasets in the manuscript following Brnich et al. (2019)  
- ACMG/AMP evidence point assignment (Tavtigian et al., 2020) based on OddsPath calibrations, ExCALIBR calibrations and gene-specific predictor calibrations 
- Scalable variant classification framework  
---

### 2. Figure Generation (`Main_Figures/`)

Contains scripts and small supporting inputs used to generate:

- **Main Figures**
- **Extended Data Figures**

Each figure directory includes:

- The plotting script  
- Required input files   
---

## Data Availability

Large supporting datasets are hosted externally on Zenodo to comply with GitHub file size limits.

**Zenodo record:** https://zenodo.org/records/18637474 

The Zenodo archive includes:

- Supplementary Data 
- Large intermediate files used for creation of the integrated variant effect dataset

Files included in this GitHub repository are:

- Analysis scripts  
- Figure-generation scripts  
- Small supporting inputs  

---

## Software Requirements

Primary languages:

- Python (3.12+)
- R (v 4.1.3)  

Additional package dependencies are specified within individual scripts.

