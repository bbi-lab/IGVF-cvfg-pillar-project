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
# 1-19, matching the "Step N" comment above each one). With no argument, all
# steps run in order followed by the flatten + condensed/expanded frame
# assembly, exactly as before. Pass a single step number as this script's
# argument to run just that one step -- normally via
# `scripts/run_variant_annotation_pipeline.sh --step N`, which also handles
# staging and env vars for you (see docs/variant_annotation_pipeline.md).
#
# prepare_gnomad_cache is a separate, optional preparation step (not part of
# the numbered step_N sequence, and never run automatically) that builds/
# refreshes Step 7's local gnomAD Hail table cache. Run it explicitly via
# `scripts/run_variant_annotation_pipeline.sh --prepare-gnomad-cache` --
# see the "gnomAD Hail table cache" prerequisite in
# docs/variant_annotation_pipeline.md for when this is (and isn't) needed.
#
# Every cvfg_variants.*/integrated_variant_effect_dataset* reference below is
# written as an explicit /work/data/... path (Dockerized steps) or
# "$VARIANT_DATA_DIR/data/..." path (the handful of plain, non-Dockerized
# shell/python commands) rather than a bare data/..., so every step reads and
# writes our own staged data/intermediate/variant_annotation/data/ regardless
# of whether the variant-annotation checkout in use happens to have its own
# files at those same relative paths -- see the "VARIANT_DATA_DIR
# path-mapping subtlety" section of docs/variant_annotation_pipeline.md.
########################################################################################################################

# Step 1: Mapping
step_1() {
src/scripts/run_map_variants.sh /work/data/cvfg_variants.0.tsv /work/data/cvfg_variants.1.tsv \
  --preferred-transcript-col preferred_transcript \
  --drop-columns target_sequence --drop-columns preferred_transcript \
  --max-clingen-concurrency 4
}

# Step 2: Replace Ensembl accessions with RefSeq
#
# --mane-file is written as /work/data/MANE.GRCh38.v1.5.summary.txt.gz rather
# than a bare data/MANE.GRCh38.v1.5.summary.txt, for the same reason as
# step_11's --requested-calibrations-file: run_remap_transcript_ids.sh remaps
# --mane-file through the same /work-vs-repo-bind-mount host-existence check
# as its positional input/output args, so a bare path would resolve against
# whichever variant-annotation checkout is in use instead of our own staged
# data/input/reference/MANE.GRCh38.v1.5.summary.txt.gz -- see
# docs/variant_annotation_pipeline.md.
step_2() {
src/scripts/run_remap_transcript_ids.sh /work/data/cvfg_variants.1.tsv /work/data/cvfg_variants.2.tsv \
  --mane-file /work/data/MANE.GRCh38.v1.5.summary.txt.gz \
  --csv-field-size-limit 10000000
}

# Step 3: Reverse translation
# Note that --wt-codon-mode unambiguous means that for WT Met and Trp "substitutions" we will generate "no_change" DNA variants in the form of codon delinses.
step_3() {
src/scripts/run_reverse_translate_protein_variants.sh /work/data/cvfg_variants.2.tsv /work/data/cvfg_variants.3.tsv \
  --include-indels \
  --wt-codon-mode unambiguous
}

# Step 4: Add VCF-style identifiers to assayed variants (both DNA and protein; already done for reverse translation candidates)
step_4() {
src/scripts/run_add_vcf_identifiers.sh /work/data/cvfg_variants.3.tsv /work/data/cvfg_variants.4.tsv --csv-field-size-limit 10000000
}

# Step 5: Add ClinGen allele IDs to reverse translations
step_5() {
src/scripts/run_add_dna_clingen_allele_ids.sh /work/data/cvfg_variants.4.tsv /work/data/cvfg_variants.5.tsv \
  --csv-field-size-limit 10000000 \
  --max-workers 5
}

# Step 6: ClinVar
step_6() {
src/scripts/run_annotate_clinvar.sh /work/data/cvfg_variants.5.tsv /work/data/cvfg_variants.6-1alt.tsv \
  --clinvar-version 201812 \
  --cache-dir ./clinvar_cache \
  --csv-field-size-limit 10000000
src/scripts/run_annotate_clinvar.sh /work/data/cvfg_variants.6-1.tsv /work/data/cvfg_variants.6-2.tsv \
  --clinvar-version 202501 \
  --cache-dir ./clinvar_cache \
  --csv-field-size-limit 10000000
src/scripts/run_annotate_clinvar.sh /work/data/cvfg_variants.6-2.tsv /work/data/cvfg_variants.6.tsv \
  --clinvar-version 202601 \
  --cache-dir ./clinvar_cache \
  --csv-field-size-limit 10000000
}

# Prepare the gnomAD Hail table cache (optional; not part of the numbered
# step_N sequence below, and never run as part of a normal pipeline run --
# see the module docstring above and the "gnomAD Hail table cache"
# prerequisite in docs/variant_annotation_pipeline.md). Run this once per
# variant-annotation checkout before Step 7's first real run against it --
# the built cache persists in that checkout's variant-annotation-gnomad-cache
# Docker volume, so an existing checkout where this has already been run
# doesn't need it again. Rebuilding takes ~6-7 hours with local[1] for a
# full gnomAD joint sites table.
prepare_gnomad_cache() {
src/scripts/run_annotate_gnomad.sh /dev/null /dev/null \
  --gnomad-version v4.1 \
  --download-only \
  --refresh-cache \
  --gnomad-ht-uri /work/gnomAD/gnomad.joint.v4.1.sites.ht
}

# Step 7: gnomAD (using local Hail table copy and Docker-volume cache; see
# prepare_gnomad_cache above for the one-time cache build/refresh)
step_7() {
# GNOMAD_CACHE_DIR=/gnomad-cache is injected automatically from the Docker volume; no --cache-dir needed.
src/scripts/run_annotate_gnomad.sh /work/data/cvfg_variants.6.tsv /work/data/cvfg_variants.7.tsv \
  --gnomad-version v4.1 \
  --require-pass \
  --callset-pass-filter any \
  --csv-field-size-limit 10000000 \
  --gnomad-ht-uri /work/gnomAD/gnomad.joint.v4.1.sites.ht \
  --log-level DEBUG
}

# Step 8: SpliceAI
step_8() {
src/scripts/run_annotate_spliceai.sh /work/data/cvfg_variants.7.tsv /work/data/cvfg_variants.8.tsv \
  --mode precomputed \
  --precomputed-snv-vcf spliceai_scores.masked.snv.hg38.vcf.gz \
  --precomputed-indel-vcf spliceai_scores.masked.indel.hg38.vcf.gz \
  --max-workers 8 \
  --csv-field-size-limit 10000000
}

# Step 9: ClinGen Evidence Repository
step_9() {
src/scripts/run_annotate_erepo.sh /work/data/cvfg_variants.8.tsv /work/data/cvfg_variants.9.tsv \
  --csv-field-size-limit 10000000
}

# Step 10: VEP mutational consequence
step_10() {
src/scripts/run_annotate_vep.sh /work/data/cvfg_variants.9.tsv /work/data/cvfg_variants.10.tsv \
  --vep-batch-size 20 \
  --row-batch-size 20 \
  --vep-timeout-seconds 60 \
  --csv-field-size-limit 10000000 \
  --log-level INFO --vep-workers 1
}

# Step 11: MaveDB variant functional classifications
#
# --requested-calibrations-file is written as /work/data/score_sets.tsv
# rather than a bare data/score_sets.tsv, for the same reason as step_13's
# predictor-file flags below: run_annotate_mavedb.sh only remaps its two
# positional input/output args through the /work-vs-repo-bind-mount
# host-existence check, so an extra flag like this one is passed through
# verbatim and would otherwise resolve against whichever variant-annotation
# checkout is in use instead of our own staged data/input/maves/score_sets.tsv
# -- see docs/variant_annotation_pipeline.md.
step_11() {
src/scripts/run_annotate_mavedb.sh /work/data/cvfg_variants.10.tsv /work/data/cvfg_variants.11.tsv \
  --requested-calibrations-file /work/data/score_sets.tsv \
  --csv-field-size-limit 10000000
}

# Step 12: Fix known MaveDB functional-classification overrides (Dockerized,
# from the CVFG pillar project rather than variant-annotation -- see
# src/scripts/run_postprocess_mavedb_functional_classifications.sh there).
step_12() {
"$CVFG_PROJECT_DIR/src/scripts/run_postprocess_mavedb_functional_classifications.sh" /work/data/cvfg_variants.11.tsv /work/data/cvfg_variants.12.tsv
}

# Step 13: REVEL and AlphaMissense
# Preparatory: regenerate the REVEL/MutPred2 training-set overlap files
# (Dockerized, from the CVFG pillar project -- see
# src/scripts/run_build_training_variant_files.sh there) before running
# annotate_predictors.
#
# All five predictor-file flags below are written as /work/data/... rather
# than bare filenames, so every input to this step comes from our own
# data/intermediate/variant_annotation/data/ (VARIANT_DATA_DIR) rather than
# whichever variant-annotation checkout happens to be in use --
# run_annotate_predictors.sh only remaps its two positional input/output
# args through the /work-vs-repo-bind-mount host-existence check; extra
# flags like these are passed through verbatim, and the container's cwd is
# the variant-annotation checkout's own bind mount (/usr/src/app), not
# /work, so a bare "data/..." path here would resolve against whichever
# checkout is in use instead of our staged data. AlphaMissense_hg38.tsv.gz
# and revel_hg38.tsv.gz (plus their .tbi indexes) must be copied into
# data/intermediate/variant_annotation/data/ before running this step --
# see docs/variant_annotation_pipeline.md.
step_13() {
"$CVFG_PROJECT_DIR/src/scripts/run_build_training_variant_files.sh"
src/scripts/run_annotate_predictors.sh /work/data/cvfg_variants.12.tsv /work/data/cvfg_variants.13.tsv \
  --alphamissense-file /work/data/AlphaMissense_hg38.tsv.gz \
  --mutpred2-properties-file /work/data/data_frame_missense_variants_MP2_properties.csv.gz \
  --revel-file /work/data/revel_hg38.tsv.gz \
  --revel-training-file /work/data/revel_training_variants.tsv \
  --mutpred2-training-file /work/data/mutpred2_training_variants.tsv \
  --csv-field-size-limit 10000000
}

# Step 14: Choose the active functional classification (Dockerized, from the
# CVFG pillar project rather than variant-annotation -- see
# src/scripts/run_add_mavedb_active_calibration_columns.sh there).
step_14() {
"$CVFG_PROJECT_DIR/src/scripts/run_add_mavedb_active_calibration_columns.sh" /work/data/cvfg_variants.13.tsv /work/data/cvfg_variants.14.tsv
}

# Step 15: Dataset names
#
# The awk pass runs directly on the host (this whole script executes with
# the variant-annotation checkout as cwd, so it isn't a bare data/... path
# here either) -- its input/output are written as
# "$VARIANT_DATA_DIR/data/..." so they land in and read from our staged
# directory regardless of cwd. The merge-columns extra-file argument is
# written as /work/data/score_sets.tsv rather than a bare data/score_sets.tsv.
# merge-columns's wrapper does remap all three of its positional arguments
# through its own /work-vs-repo-bind-mount host-existence check, but that
# check would still incorrectly prefer a variant-annotation checkout's own
# data/score_sets.tsv over our staged data/input/maves/score_sets.tsv if one
# happens to exist there -- so it's forced to /work explicitly, same as
# step_11 and step_13 -- see docs/variant_annotation_pipeline.md.
step_15() {
awk -F'\t' -v OFS='\t' '
  NR==1 { sub(/\r$/, ""); for (i=1; i<=NF; i++) if ($i=="variant_urn") col=i; print $0, "score_set_urn" }
  NR>1  { sub(/\r$/, ""); val = col ? $col : ""; sub(/#.*$/, "", val); print $0, val }
' "$VARIANT_DATA_DIR/data/cvfg_variants.14.tsv" > "$VARIANT_DATA_DIR/data/cvfg_variants.15.temp.tsv"
src/scripts/run_utilities.sh merge-columns \
  /work/data/cvfg_variants.15.temp.tsv \
  /work/data/score_sets.tsv \
  /work/data/cvfg_variants.15.tsv \
  --key-col "score_set_urn:score_set_urn" \
  --add-col "dataset_name" \
  --csv-field-size-limit 10000000
rm "$VARIANT_DATA_DIR/data/cvfg_variants.15.temp.tsv"
}

# Step 16: Assay metadata from Supplementary Data 3
#
# The extra-file argument is written as /work/data/Supplementary_Data_3.xlsx
# rather than a bare data/Supplementary_Data_3.xlsx, for the same reason as
# step_15's score_sets.tsv above: Supplementary_Data_3.xlsx lives in this
# project's own data/input/maves/ (not a variant-annotation checkout), and
# merge-columns's wrapper's host-existence check on a bare path would still
# incorrectly prefer a variant-annotation checkout's own copy if one happens
# to exist there -- see docs/variant_annotation_pipeline.md.
step_16() {
src/scripts/run_utilities.sh merge-columns \
  /work/data/cvfg_variants.15.tsv \
  /work/data/Supplementary_Data_3.xlsx \
  /work/data/cvfg_variants.16.tsv \
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

# Step 17: Simplified consequence (Dockerized, from the CVFG pillar project
# rather than variant-annotation -- see
# src/scripts/run_annotate_simplified_consequence.sh there).
step_17() {
"$CVFG_PROJECT_DIR/src/scripts/run_annotate_simplified_consequence.sh" /work/data/cvfg_variants.16.tsv /work/data/cvfg_variants.17.tsv
}

# Step 18: Recalculate ClinGen classification without functional-assay evidence
# (Dockerized, from the CVFG pillar project -- see
# src/scripts/run_recalculate_clingen_classification.sh there).
step_18() {
"$CVFG_PROJECT_DIR/src/scripts/run_recalculate_clingen_classification.sh" /work/data/cvfg_variants.17.tsv /work/data/cvfg_variants.18.tsv
}

# Step 19: Flag variants (Dockerized, from the CVFG pillar project rather
# than variant-annotation -- see src/scripts/run_flag_variants.sh there).
step_19() {
"$CVFG_PROJECT_DIR/src/scripts/run_flag_variants.sh" /work/data/cvfg_variants.18.tsv /work/data/cvfg_variants.19.tsv
}

LAST_STEP=19

# Flatten, then assemble the condensed and expanded data frames. Not part of
# the numbered step_N sequence (see module docstring above) -- always runs in
# full, as the tail of a complete pipeline run.
build_condensed_and_expanded_frames() {

# Flatten
src/scripts/run_flatten_dna_variants.sh /work/data/cvfg_variants.19.tsv /work/data/cvfg_variants.19.flat.tsv

########################################################################################################################
# Condensed data frame
########################################################################################################################

# Filter out unmapped variants.
src/scripts/run_utilities.sh filter-rows \
  /work/data/cvfg_variants.19.tsv \
  /work/data/cvfg_variants.19.mapped.tsv \
  --value-col "mapped_hgvs_g,mapped_hgvs_c,mapped_hgvs_p" \
  --match any \
  --csv-field-size-limit 10000000

# Filter, reorder, and rename columns to match the CVFG dataframe format.
src/scripts/run_utilities.sh rename-columns \
  /work/data/cvfg_variants.19.mapped.tsv \
  /work/data/integrated_variant_effect_dataset.condensed.tsv \
  --csv-field-size-limit 10000000 \
  --reorder \
  --keep-col "dataset_name:Dataset" \
  --keep-col "gene_symbol:Gene" \
  --keep-col "gene_hgnc_id:HGNC_ID" \
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
  --keep-col "clingen_evidence_repository.Disease:Disease_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Mondo Id:Mondo Id_ClinGen_repo" \
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
' "$VARIANT_DATA_DIR/data/integrated_variant_effect_dataset.condensed.tsv" > "$VARIANT_DATA_DIR/data/file.tsv" && mv "$VARIANT_DATA_DIR/data/file.tsv" "$VARIANT_DATA_DIR/data/integrated_variant_effect_dataset.condensed.tsv"

########################################################################################################################
# Expanded data frame
########################################################################################################################

src/scripts/run_utilities.sh filter-rows \
  /work/data/cvfg_variants.19.flat.tsv \
  /work/data/cvfg_variants.19.flat.mapped.tsv \
  --value-col "mapped_hgvs_g,mapped_hgvs_c,mapped_hgvs_p" \
  --match any \
  --csv-field-size-limit 10000000

# Filter, reorder, and rename columns to match the CVFG dataframe format.
src/scripts/run_utilities.sh rename-columns \
  /work/data/cvfg_variants.19.flat.mapped.tsv \
  /work/data/integrated_variant_effect_dataset.tsv \
  --csv-field-size-limit 10000000 \
  --reorder \
  --keep-col "dataset_name:Dataset" \
  --keep-col "gene_symbol:Gene" \
  --keep-col "gene_hgnc_id:HGNC_ID" \
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
  --keep-col "clingen_evidence_repository.Disease:Disease_ClinGen_repo" \
  --keep-col "clingen_evidence_repository.Mondo Id:Mondo Id_ClinGen_repo" \
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
' "$VARIANT_DATA_DIR/data/integrated_variant_effect_dataset.tsv" > "$VARIANT_DATA_DIR/data/file.tsv" && mv "$VARIANT_DATA_DIR/data/file.tsv" "$VARIANT_DATA_DIR/data/integrated_variant_effect_dataset.tsv"

}

########################################################################################################################
# Entry point: run everything, just one requested step, or the optional
# prepare-gnomad-cache preparation step.
########################################################################################################################

requested_step="${1:-}"

if [[ "$requested_step" == "prepare-gnomad-cache" ]]; then
  prepare_gnomad_cache
  exit 0
fi

if [[ -n "$requested_step" ]]; then
  step_fn="step_${requested_step}"
  if ! declare -F "$step_fn" > /dev/null; then
    echo "error: no such step '$requested_step' (valid: 1-$LAST_STEP, or 'prepare-gnomad-cache')" >&2
    exit 1
  fi
  "$step_fn"
else
  for ((n = 1; n <= LAST_STEP; n++)); do
    "step_$n"
  done
  build_condensed_and_expanded_frames
fi
