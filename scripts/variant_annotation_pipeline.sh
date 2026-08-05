########################################################################################################################
# IGVF CVFG pipeline
#
# Not meant to be run directly. Invoke via
# scripts/run_variant_annotation_pipeline.sh, which:
#   - stages data/raw_mave_data/ into the shared data dir this script reads
#     data/... paths from,
#   - cd's into the variant-annotation checkout (so the src/scripts/run_*.sh
#     wrapper calls below resolve) and exports VARIANT_DATA_DIR to point at
#     that staged data dir,
#   - exports CVFG_PROJECT_DIR (this repo's root) so the "Flag variants" step
#     below can find its own Dockerized wrapper regardless of cwd,
#   - and gzips the two final integrated_variant_effect_dataset files into
#     data/mave_data/ once this script finishes.
# See docs/variant_annotation_pipeline.md for the full data flow and why
# VARIANT_DATA_DIR is repurposed this way.
#
# Steps are defined as step_N functions below: step_N reads
# data/cvfg_variants.<N-1>.tsv and writes data/cvfg_variants.<N>.tsv (N is
# 1-18, matching the "Step N" comment above each one). With no argument, all
# steps run in order followed by the flatten + condensed/expanded frame
# assembly, exactly as before. Pass a single step number as this script's
# argument to run just that one step -- normally via
# `scripts/run_variant_annotation_pipeline.sh --step N`, which also handles
# staging and env vars for you (see docs/variant_annotation_pipeline.md).
########################################################################################################################

# Step 1: Mapping
step_1() {
src/scripts/run_map_variants.sh data/cvfg_variants.0.tsv data/cvfg_variants.1.tsv \
  --preferred-transcript-col preferred_transcript \
  --drop-columns target_sequence --drop-columns preferred_transcript \
  --max-clingen-concurrency 4
}

# Step 2: Replace Ensembl accessions with RefSeq
step_2() {
src/scripts/run_remap_transcript_ids.sh data/cvfg_variants.1.tsv data/cvfg_variants.2.tsv \
  --mane-file data/MANE.GRCh38.v1.5.summary.txt \
  --csv-field-size-limit 10000000
}

# Step 3: Reverse translation
# Note that --wt-codon-mode unambiguous means that for WT Met and Trp "substitutions" we will generate "no_change" DNA variants in the form of codon delinses.
step_3() {
src/scripts/run_reverse_translate_protein_variants.sh data/cvfg_variants.2.tsv data/cvfg_variants.3.tsv \
  --include-indels \
  --wt-codon-mode unambiguous
}

# Step 4: Add VCF-style identifiers to assayed variants (both DNA and protein; already done for reverse translation candidates)
step_4() {
src/scripts/run_add_vcf_identifiers.sh data/cvfg_variants.3.tsv data/cvfg_variants.4.tsv --csv-field-size-limit 10000000
}

# Step 5: Add ClinGen allele IDs to reverse translations
step_5() {
src/scripts/run_add_dna_clingen_allele_ids.sh data/cvfg_variants.4.tsv data/cvfg_variants.5.tsv \
  --csv-field-size-limit 10000000 \
  --max-workers 5
}

# Step 6: ClinVar
step_6() {
src/scripts/run_annotate_clinvar.sh data/cvfg_variants.5.tsv data/cvfg_variants.6-1alt.tsv \
  --clinvar-version 201812 \
  --cache-dir ./clinvar_cache \
  --csv-field-size-limit 10000000
src/scripts/run_annotate_clinvar.sh data/cvfg_variants.6-1.tsv data/cvfg_variants.6-2.tsv \
  --clinvar-version 202501 \
  --cache-dir ./clinvar_cache \
  --csv-field-size-limit 10000000
src/scripts/run_annotate_clinvar.sh data/cvfg_variants.6-2.tsv data/cvfg_variants.6.tsv \
  --clinvar-version 202601 \
  --cache-dir ./clinvar_cache \
  --csv-field-size-limit 10000000
}

# Step 7: gnomAD (using local Hail table copy and Docker-volume cache)
# Build the cache first if needed (one-time; ~6-7 hours with local[1]):
step_7() {
src/scripts/run_annotate_gnomad.sh /dev/null /dev/null \
  --gnomad-version v4.1 \
  --download-only \
  --refresh-cache \
  --gnomad-ht-uri /work/gnomAD/gnomad.joint.v4.1.sites.ht
# GNOMAD_CACHE_DIR=/gnomad-cache is injected automatically from the Docker volume; no --cache-dir needed.
src/scripts/run_annotate_gnomad.sh data/cvfg_variants.6.tsv data/cvfg_variants.7.tsv \
  --gnomad-version v4.1 \
  --require-pass \
  --callset-pass-filter any \
  --csv-field-size-limit 10000000 \
  --gnomad-ht-uri /work/gnomAD/gnomad.joint.v4.1.sites.ht \
  --log-level DEBUG
}

# Step 8: SpliceAI
step_8() {
src/scripts/run_annotate_spliceai.sh data/cvfg_variants.7.tsv data/cvfg_variants.8.tsv \
  --mode precomputed \
  --precomputed-snv-vcf spliceai_scores.masked.snv.hg38.vcf.gz \
  --precomputed-indel-vcf spliceai_scores.masked.indel.hg38.vcf.gz \
  --max-workers 8 \
  --csv-field-size-limit 10000000
}

# Step 9: ClinGen Evidence Repository
step_9() {
src/scripts/run_annotate_erepo.sh data/cvfg_variants.8.tsv data/cvfg_variants.9.tsv \
  --csv-field-size-limit 10000000
}

# Step 10: VEP mutational consequence
step_10() {
src/scripts/run_annotate_vep.sh data/cvfg_variants.9.tsv data/cvfg_variants.10.tsv \
  --vep-batch-size 20 \
  --row-batch-size 20 \
  --vep-timeout-seconds 60 \
  --csv-field-size-limit 10000000 \
  --log-level INFO --vep-workers 1
}

# Step 11: MaveDB variant functional classifications
step_11() {
src/scripts/run_annotate_mavedb.sh data/cvfg_variants.10.tsv data/cvfg_variants.11.tsv \
  --requested-calibrations-file data/score_sets.tsv \
  --csv-field-size-limit 10000000
}

# Step 12: REVEL and AlphaMissense
# Preparatory: regenerate the REVEL/MutPred2 training-set overlap files
# (Dockerized, from the CVFG pillar project -- see
# src/scripts/run_build_training_variant_files.sh there) before running
# annotate_predictors.
step_12() {
"$CVFG_PROJECT_DIR/src/scripts/run_build_training_variant_files.sh"
src/scripts/run_annotate_predictors.sh data/cvfg_variants.11.tsv data/cvfg_variants.12.tsv \
  --alphamissense-file AlphaMissense_hg38.tsv.gz \
  --mutpred2-properties-file data/data_frame_missense_variants_MP2_properties.csv.gz \
  --revel-file revel_hg38.tsv.gz \
  --revel-training-file data/revel_training_variants.tsv \
  --mutpred2-training-file data/mutpred2_training_variants.tsv \
  --csv-field-size-limit 10000000
}

# Step 13: Choose the active functional classification
step_13() {
src/scripts/add_mavedb_active_calibration_columns.sh
}

# Step 14: Dataset names
step_14() {
awk -F'\t' -v OFS='\t' '
  NR==1 { sub(/\r$/, ""); for (i=1; i<=NF; i++) if ($i=="variant_urn") col=i; print $0, "score_set_urn" }
  NR>1  { sub(/\r$/, ""); val = col ? $col : ""; sub(/#.*$/, "", val); print $0, val }
' data/cvfg_variants.13.tsv > data/cvfg_variants.14.temp.tsv
src/scripts/run_utilities.sh merge-columns \
  data/cvfg_variants.14.temp.tsv \
  data/score_sets.tsv \
  data/cvfg_variants.14.tsv \
  --key-col "score_set_urn:score_set_urn" \
  --add-col "dataset_name" \
  --csv-field-size-limit 10000000
rm data/cvfg_variants.14.temp.tsv
}

# Step 15: Assay metadata from Supplementary Data 3
step_15() {
src/scripts/run_utilities.sh merge-columns \
  data/cvfg_variants.14.tsv \
  data/Supplementary_Data_3.xlsx \
  data/cvfg_variants.15.tsv \
  --csv-field-size-limit 10000000 \
  --extra-worksheet Curation \
  --key-col "dataset_name:Dataset Name" \
  --add-col "Gene:gene_symbol" \
  --add-col "HGNC ID:gene_hgnc_id" \
  --add-col "Detects Splicing Variants?:assay_detects_splicing_effects" \
  --add-col "Author Provided Transcript ID:author_provided_transcript_id" \
  --add-col "Ensembl Transcript ID:ensembl_transcript_id_from_assay_metadata" \
  --add-col "RefSeq Transcript ID:refseq_transcript_id_from_assay_metadata" \
  --add-col "Interval 1 Name" \
  --add-col "Interval 1 Range" \
  --add-col "Interval 1 Class" \
  --add-col "Interval 2 Name" \
  --add-col "Interval 2 Range" \
  --add-col "Interval 2 Class" \
  --add-col "Interval 3 Name" \
  --add-col "Interval 3 Range" \
  --add-col "Interval 3 Class" \
  --add-col "Interval 4 Name" \
  --add-col "Interval 4 Range" \
  --add-col "Interval 4 Class" \
  --add-col "Interval 5 Name" \
  --add-col "Interval 5 Range" \
  --add-col "Interval 5 Class" \
  --add-col "Interval 6 Name" \
  --add-col "Interval 6 Range" \
  --add-col "Interval 6 Class"
}

# Step 16: Extra fields for GMM
step_16() {
python3 -m src.annotate_simplified_consequence data/cvfg_variants.15.tsv data/cvfg_variants.16.tsv \
  --consequence-map-file data/extended_ensembl_consequence.csv.gz \
  --csv-field-size-limit 10000000
}

# Step 17: Recalculate ClinGen classification without functional-assay evidence
# (Dockerized, from the CVFG pillar project -- see
# src/scripts/run_recalculate_clingen_classification.sh there).
step_17() {
"$CVFG_PROJECT_DIR/src/scripts/run_recalculate_clingen_classification.sh" data/cvfg_variants.16.tsv data/cvfg_variants.17.tsv
}

# Step 18: Flag variants (Dockerized, from the CVFG pillar project rather
# than variant-annotation -- see src/scripts/run_flag_variants.sh there).
step_18() {
"$CVFG_PROJECT_DIR/src/scripts/run_flag_variants.sh" data/cvfg_variants.17.tsv data/cvfg_variants.18.tsv
}

LAST_STEP=18

# Flatten, then assemble the condensed and expanded data frames. Not part of
# the numbered step_N sequence (see module docstring above) -- always runs in
# full, as the tail of a complete pipeline run.
build_condensed_and_expanded_frames() {

# Flatten
src/scripts/run_flatten_dna_variants.sh data/cvfg_variants.18.tsv data/cvfg_variants.18.flat.tsv

########################################################################################################################
# Condensed data frame
########################################################################################################################

# Filter out unmapped variants.
src/scripts/run_utilities.sh filter-rows \
  data/cvfg_variants.18.tsv \
  data/cvfg_variants.18.mapped.tsv \
  --value-col "mapped_hgvs_g,mapped_hgvs_c,mapped_hgvs_p" \
  --match any \
  --csv-field-size-limit 10000000

# Filter, reorder, and rename columns to match the CVFG dataframe format.
src/scripts/run_utilities.sh rename-columns \
  data/cvfg_variants.18.mapped.tsv \
  data/integrated_variant_effect_dataset.condensed.tsv \
  --csv-field-size-limit 10000000 \
  --reorder \
  --keep-col "dataset_name:Dataset" \
  --keep-col "gene_symbol:Gene" \
  --keep-col "variant_urn:mavedb_variant_urn" \
  --keep-col "mapped_hgvs_g_chromosome:Chrom" \
  --keep-col "strand:Strand" \
  --keep-col "mapped_hgvs_g_start:hg38_start" \
  --keep-col "mapped_hgvs_g_stop:hg38_end" \
  --keep-col "mapped_hgvs_g_ref:ref_allele" \
  --keep-col "mapped_hgvs_g_alt:alt_allele" \
  --keep-col "author_provided_transcript_id:auth_transcript_id" \
  --keep-col "mapped_hgvs_c_start:transcript_pos" \
  --keep-col "mapped_hgvs_c_ref:transcript_ref" \
  --keep-col "mapped_hgvs_c_alt:transcript_alt" \
  --keep-col "mapped_hgvs_p_start:aa_pos" \
  --keep-col "mapped_hgvs_p_ref:aa_ref" \
  --keep-col "mapped_hgvs_p_alt:aa_alt" \
  --keep-col "mapped_hgvs_c:hgvs_c" \
  --keep-col "mapped_hgvs_p:hgvs_p" \
  --keep-col "vep.mutational_consequences:consequence" \
  --keep-col "vep.most_severe_mutational_consequence:most_severe_mutational_consequence" \
  --keep-col "score:auth_reported_score" \
  --keep-col "rna_score" \
  --keep-col "rna_score_d6" \
  --keep-col "rna_score_d20" \
  --keep-col "mavedb.active_calibration.functional_class_label:auth_reported_func_class" \
  --keep-col "mavedb.active_calibration.functional_classification:auth_reported_func_class_category" \
  --keep-col "assay_detects_splicing_effects:splice_measure" \
  --keep-col "gnomad.v4_1.minor_allele_frequency:gnomad_MAF" \
  --keep-col "clinvar.202601.clinical_significance:clinvar_sig_2026" \
  --keep-col "clinvar.202601.review_status:clinvar_star_2026" \
  --keep-col "clinvar.202601.last_review_date:clinvar_date_last_reviewed_2026" \
  --keep-col "clinvar.202501.clinical_significance:clinvar_sig_2025" \
  --keep-col "clinvar.202501.review_status:clinvar_star_2025" \
  --keep-col "clinvar.202501.last_review_date:clinvar_date_last_reviewed_2025" \
  --keep-col "clinvar.201812.clinical_significance:clinvar_sig_2018" \
  --keep-col "clinvar.201812.review_status:clinvar_star_2018" \
  --keep-col "clinvar.201812.last_review_date:clinvar_date_last_reviewed_2018" \
  --keep-col "assayed_variant_level:nucleotide_or_aa" \
  --keep-col "ensembl_transcript_id_from_assay_metadata:Ensembl Transcript ID" \
  --keep-col "refseq_transcript_id_from_assay_metadata:RefSeq Transcript ID" \
  --keep-col "Interval 1 Name" \
  --keep-col "Interval 1 Range" \
  --keep-col "Interval 1 Class" \
  --keep-col "Interval 2 Name" \
  --keep-col "Interval 2 Range" \
  --keep-col "Interval 2 Class" \
  --keep-col "Interval 3 Name" \
  --keep-col "Interval 3 Range" \
  --keep-col "Interval 3 Class" \
  --keep-col "Interval 4 Name" \
  --keep-col "Interval 4 Range" \
  --keep-col "Interval 4 Class" \
  --keep-col "Interval 5 Name" \
  --keep-col "Interval 5 Range" \
  --keep-col "Interval 5 Class" \
  --keep-col "Interval 6 Name" \
  --keep-col "Interval 6 Range" \
  --keep-col "Interval 6 Class" \
  --keep-col "spliceai.ds_ag:spliceAI_DS_AG" \
  --keep-col "spliceai.ds_al:spliceAI_DS_AL" \
  --keep-col "spliceai.ds_dg:spliceAI_DS_DG" \
  --keep-col "spliceai.ds_dl:spliceAI_DS_DL" \
  --keep-col "spliceai.dp_ag:spliceAI_DP_AG" \
  --keep-col "spliceai.dp_al:spliceAI_DP_AL" \
  --keep-col "spliceai.dp_dg:spliceAI_DP_DG" \
  --keep-col "spliceai.dp_dl:spliceAI_DP_DL" \
  --keep-col "spliceai.dp_dl:spliceAI_DP_DL" \
  --keep-col "clingen_evidence_repository.ClinVar Variation Id:ClinVar Variation Id_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Allele Registry Id:Allele Registry Id_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Disease Mondo Id:Mondo Id_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Mode of Inheritance:Mode of Inheritance_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Assertion:Assertion_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Applied Evidence Codes (Met):Applied Evidence Codes (Met)_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Applied Evidence Codes (Not Met):Applied Evidence Codes (Not Met)_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Summary of interpretation:Summary of interpretation_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.PubMed Articles:PubMed Articles_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Expert Panel:Expert Panel_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Guideline:Guideline_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Approval Date:Approval Date_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Published Date:Published Date_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Retracted:Retracted_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Evidence Repo Link:Evidence Repo Link_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Uuid:Uuid_ClinGen_repo" \
  --keep-col "Updated_Classification_ClinGen_repo" \
  --keep-col "Updated_Evidence Codes_ClinGen_repo" \
  --keep-col "revel.score:REVEL" \
  --keep-col "alphamissense.pathogenicity:AM_score" \
  --keep-col "alphamissense.class:AM_class" \
  --keep-col "mutpred2.score:MutPred2" \
  --keep-col "simplified_consequence" \
  --keep-col "condensed_consequence" \
  --keep-col "splice_variant" \
  --keep-col "splice_var_amino" \
  --keep-col "Flag"

# Change codes used in the nucleotide_or_aa column.
awk -F'\t' -v OFS='\t' '
  NR==1 { for (i=1; i<=NF; i++) if ($i=="nucleotide_or_aa") col=i; print }
  NR>1  { if (col) { if ($col=="protein") $col="aa"; else if ($col=="dna") $col="nt" }; print }
' data/integrated_variant_effect_dataset.condensed.tsv > data/file.tsv && mv data/file.tsv data/integrated_variant_effect_dataset.condensed.tsv

# Fix two things
src/scripts/postprocess_integrated_variant_effect_dataset.sh data/integrated_variant_effect_dataset.condensed.tsv data/file.tsv && mv data/file.tsv data/integrated_variant_effect_dataset.condensed.tsv

########################################################################################################################
# Expanded data frame
########################################################################################################################

src/scripts/run_utilities.sh filter-rows \
  data/cvfg_variants.18.flat.tsv \
  data/cvfg_variants.18.flat.mapped.tsv \
  --value-col "mapped_hgvs_g,mapped_hgvs_c,mapped_hgvs_p" \
  --match any \
  --csv-field-size-limit 10000000

# Filter, reorder, and rename columns to match the CVFG dataframe format.
src/scripts/run_utilities.sh rename-columns \
  data/cvfg_variants.18.flat.mapped.tsv \
  data/integrated_variant_effect_dataset.tsv \
  --csv-field-size-limit 10000000 \
  --reorder \
  --keep-col "dataset_name:Dataset" \
  --keep-col "gene_symbol:Gene" \
  --keep-col "variant_urn:mavedb_variant_urn" \
  --keep-col "mapped_hgvs_g_chromosome:Chrom" \
  --keep-col "strand:Strand" \
  --keep-col "mapped_hgvs_g_start:hg38_start" \
  --keep-col "mapped_hgvs_g_stop:hg38_end" \
  --keep-col "mapped_hgvs_g_ref:ref_allele" \
  --keep-col "mapped_hgvs_g_alt:alt_allele" \
  --keep-col "author_provided_transcript_id:auth_transcript_id" \
  --keep-col "mapped_hgvs_c_start:transcript_pos" \
  --keep-col "mapped_hgvs_c_ref:transcript_ref" \
  --keep-col "mapped_hgvs_c_alt:transcript_alt" \
  --keep-col "mapped_hgvs_p_start:aa_pos" \
  --keep-col "mapped_hgvs_p_ref:aa_ref" \
  --keep-col "mapped_hgvs_p_alt:aa_alt" \
  --keep-col "mapped_hgvs_c:hgvs_c" \
  --keep-col "mapped_hgvs_p:hgvs_p" \
  --keep-col "vep.mutational_consequences:consequence" \
  --keep-col "vep.most_severe_mutational_consequence:most_severe_mutational_consequence" \
  --keep-col "score:auth_reported_score" \
  --keep-col "rna_score" \
  --keep-col "rna_score_d6" \
  --keep-col "rna_score_d20" \
  --keep-col "mavedb.active_calibration.functional_class_label:auth_reported_func_class" \
  --keep-col "mavedb.active_calibration.functional_classification:auth_reported_func_class_category" \
  --keep-col "assay_detects_splicing_effects:splice_measure" \
  --keep-col "gnomad.v4_1.minor_allele_frequency:gnomad_MAF" \
  --keep-col "clinvar.202601.clinical_significance:clinvar_sig_2026" \
  --keep-col "clinvar.202601.review_status:clinvar_star_2026" \
  --keep-col "clinvar.202601.last_review_date:clinvar_date_last_reviewed_2026" \
  --keep-col "clinvar.202501.clinical_significance:clinvar_sig_2025" \
  --keep-col "clinvar.202501.review_status:clinvar_star_2025" \
  --keep-col "clinvar.202501.last_review_date:clinvar_date_last_reviewed_2025" \
  --keep-col "clinvar.201812.clinical_significance:clinvar_sig_2018" \
  --keep-col "clinvar.201812.review_status:clinvar_star_2018" \
  --keep-col "clinvar.201812.last_review_date:clinvar_date_last_reviewed_2018" \
  --keep-col "assayed_variant_level:nucleotide_or_aa" \
  --keep-col "ensembl_transcript_id_from_assay_metadata:Ensembl Transcript ID" \
  --keep-col "refseq_transcript_id_from_assay_metadata:RefSeq Transcript ID" \
  --keep-col "Interval 1 Name" \
  --keep-col "Interval 1 Range" \
  --keep-col "Interval 1 Class" \
  --keep-col "Interval 2 Name" \
  --keep-col "Interval 2 Range" \
  --keep-col "Interval 2 Class" \
  --keep-col "Interval 3 Name" \
  --keep-col "Interval 3 Range" \
  --keep-col "Interval 3 Class" \
  --keep-col "Interval 4 Name" \
  --keep-col "Interval 4 Range" \
  --keep-col "Interval 4 Class" \
  --keep-col "Interval 5 Name" \
  --keep-col "Interval 5 Range" \
  --keep-col "Interval 5 Class" \
  --keep-col "Interval 6 Name" \
  --keep-col "Interval 6 Range" \
  --keep-col "Interval 6 Class" \
  --keep-col "spliceai.ds_ag:spliceAI_DS_AG" \
  --keep-col "spliceai.ds_al:spliceAI_DS_AL" \
  --keep-col "spliceai.ds_dg:spliceAI_DS_DG" \
  --keep-col "spliceai.ds_dl:spliceAI_DS_DL" \
  --keep-col "spliceai.dp_ag:spliceAI_DP_AG" \
  --keep-col "spliceai.dp_al:spliceAI_DP_AL" \
  --keep-col "spliceai.dp_dg:spliceAI_DP_DG" \
  --keep-col "spliceai.dp_dl:spliceAI_DP_DL" \
  --keep-col "spliceai.dp_dl:spliceAI_DP_DL" \
  --keep-col "clingen_evidence_repository.ClinVar Variation Id:ClinVar Variation Id_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Allele Registry Id:Allele Registry Id_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Disease Mondo Id:Mondo Id_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Mode of Inheritance:Mode of Inheritance_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Assertion:Assertion_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Applied Evidence Codes (Met):Applied Evidence Codes (Met)_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Applied Evidence Codes (Not Met):Applied Evidence Codes (Not Met)_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Summary of interpretation:Summary of interpretation_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.PubMed Articles:PubMed Articles_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Expert Panel:Expert Panel_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Guideline:Guideline_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Approval Date:Approval Date_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Published Date:Published Date_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Retracted:Retracted_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Evidence Repo Link:Evidence Repo Link_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Uuid:Uuid_ClinGen_repo" \
  --keep-col "Updated_Classification_ClinGen_repo" \
  --keep-col "Updated_Evidence Codes_ClinGen_repo" \
  --keep-col "revel.score:REVEL" \
  --keep-col "alphamissense.pathogenicity:AM_score" \
  --keep-col "alphamissense.class:AM_class" \
  --keep-col "mutpred2.score:MutPred2" \
  --keep-col "simplified_consequence" \
  --keep-col "condensed_consequence" \
  --keep-col "splice_variant" \
  --keep-col "splice_var_amino" \
  --keep-col "Flag"

# Change codes used in the nucleotide_or_aa column.
awk -F'\t' -v OFS='\t' '
  NR==1 { for (i=1; i<=NF; i++) if ($i=="nucleotide_or_aa") col=i; print }
  NR>1  { if (col) { if ($col=="protein") $col="aa"; else if ($col=="dna") $col="nt" }; print }
' data/integrated_variant_effect_dataset.tsv > data/file.tsv && mv data/file.tsv data/integrated_variant_effect_dataset.tsv

# Fix two things
src/scripts/postprocess_integrated_variant_effect_dataset.sh data/integrated_variant_effect_dataset.tsv data/file.tsv && mv data/file.tsv data/integrated_variant_effect_dataset.tsv

}

########################################################################################################################
# Entry point: run everything, or just one requested step.
########################################################################################################################

requested_step="${1:-}"

if [[ -n "$requested_step" ]]; then
  step_fn="step_${requested_step}"
  if ! declare -F "$step_fn" > /dev/null; then
    echo "error: no such step '$requested_step' (valid: 1-$LAST_STEP)" >&2
    exit 1
  fi
  "$step_fn"
else
  for ((n = 1; n <= LAST_STEP; n++)); do
    "step_$n"
  done
  build_condensed_and_expanded_frames
fi
