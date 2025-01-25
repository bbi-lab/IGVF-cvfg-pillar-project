# Pillar Project dataframe columns definitions

## This file provides definitions for all the columns in the pillar project dataframe and other important information. Please read carefully before using the dataframe for any analysis.

## **Important** This DataFrame contains score sets that are evaluated at *both* the ***nucleotide*** and ***amino acid levels***. Amino acid-level score sets are mapped back to nucleotide-level information using the VEP web interface. Since each amino acid-level variant can correspond to multiple nucleotide substitutions, the DataFrame is structured to maintain these relationships while preserving the original scores provided by the authors. To achieve this, the DataFrame adopts a "short format," where each amino acid-level variant’s corresponding nucleotide-level information is concatenated and separated by the '^' character. By maintaining this structure, the DataFrame avoids confusion and upholds the integrity of the authors' original data. If parsing or separating the nucleotide-level information is needed, the '^' character can be used as a delimiter. 


### Clinvar version: October 15th, 2024 

### gnomad v4

#### Column definitions

| **Column name**           | **Datatype** | **Description**
|---------------------------|--------------|-------------------------

| ID                        | object       | Unique identifier given to each author-provided variant in the dataframe

| Dataset                   | object       | Name given to each individual scoreset, will match MAVE curation sheet

| Gene                      | object       | HGNC gene symbol

| HGNC_id                   | int64        | HGNC gene ID

| Chrom                     | object       | Chromosome 

| hg19_pos                  | float64      | hg19 (GRCh37) assembly genomic coordinates

| hg38_start                | float64      | hg38 (GRCh38) assembly genomic coordinates, start position of variant

| hg38_end                  | float64      | hg38 (GRCh38) assembly genomic coordinates, end position of variant

| ref_allele                | object       | reference allele 

| alt_allele                | object       | alternate allele

|auth_transcript_id         | object       | author provided transcript ID

|transcript_pos             | object       | transcript position

|transcript_ref             | object       | transcript reference allele 

|transcript_alt             | object       | transcript alternate allele

|aa_pos                     | object       | amino acid position

|aa_ref                     | object       | amino acid reference protein

|aa_alt                     | object       | amino acid alternate protein

|hgvs_c                     | object       | hgvs c. nomenclature

|hgvs_p                     | object       | hgvs p. nomenclature

|consequence                | object       | VEP annotate consequence

|auth_reported_score        | float64      | author reported score that is used for analyses in their paper; could be normalized or statistically manipulated (not necessarily the raw score from the experiment or the average of replicates)

|auth_reported_rep_score    | object       | replicate scores provided by author (could be raw scores) separated by ';'

|auth_reported_func_class   | object       | functional class of the variant reported by the authors

|auth_reported_func_class   | object       | functional class of the variant reported by the authors

|auth_reported_normal_min   | float64      | normal functional class minimum score

|auth_reported_normal_max   | float64      | normal functional class maximum score

|auth_reported_abnormal_min | float64      | abnormal functional class minimum score

|auth_reported_abnormal_max | float64      | abnormal functional class maximum score

|splice_measure             | object       | does the assay measure splicing; 'Yes' or 'No' values

|gnomad_MAF                 | float64      | gnomad v4 MAF

|clinvar_sig                | object       | Clinvar siginificance; last download October 15, 2024

|clinvar_star               | float64      | if found in clinvar how many stars does this variant get

|clinvar_date_last_reviewed | float64      | if found in clinvar, when was the variant last reviewed

|nucleotide_or_aa           | object       | is this assay performed on the nucleotide level or the amino acid level? values are 'nucleotide' or 'aa'

|MaveDB URN                 | object       | MaveDB URN number (if this score set is in MaveDB)

|Ensembl_transcript_ID      | object       | Matches author provided transcript ID, if none provided, the assay scoreset is matched to a transcript ID

|Ref_seq_transcript_ID      | object       | Matches author provided transcript ID, if none provided, the assay scoreset is matched to a transcript ID

|Model_system               | object       | Model system the assay was performed in

|Assay_type                 | object       | One of cell_viability, reporter, or direct protein function

|Phenotype_measured         | object       | Phenotype measured in the assay

|Phenotype_detail           | object       | A more detailed description of the phenotype measured in the assay

|IGVF_produced              | object       | Produced by IGVF? Values are 'Yes' or 'No'

|Flag                       | object       | Variants marked with `*` are flagged as having incorrect VEP annotations for consequence 
