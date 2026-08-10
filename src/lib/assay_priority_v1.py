"""Deduplication priority order for MAVE assays that share a gene.

`ASSAY_PRIORITY_LIST` is used by the OddsPath calibration pipeline
(`notebooks/analysis/OddsPath_classifications.ipynb` and
`notebooks/analysis/Variant_Classification_analysis.ipynb`) to resolve variants that
have functional scores from more than one assay covering the same gene.
Each call site builds `{name: i for i, name in enumerate(ASSAY_PRIORITY_LIST)}`,
maps it onto a `Dataset` column (assays not in the list fall back to
`9999`), then keeps only the lowest-priority-index row per variant via
`sort_values("assay_priority").drop_duplicates(..., keep="first")`.

This list matches the priority order used in the original preprint
submission (the `"v1"` deduplication strategy -- see
`docs/variant_classification.md`) and nothing else. It previously also
carried two further sections of unreviewed, best-guess orderings for
assays added before and after the preprint submission; those were removed
because they represented guesses that were never actually reviewed or
decided, not a real methodology choice, and their presence made `"v1"`
inconsistent with what the preprint submission actually did. Every assay
not in this list -- which today includes all `BRCA2` assays other than
`BRCA2_Hu_2024`, all `LDLR`, `CBS`, `HMBS`, and `PAX6` assays, and several
others -- now falls back to the same unranked `9999` priority as any other
unlisted assay, so ties among them are resolved by an unstable sort rather
than any considered order. See
`docs/assay_priority_questions.md` (section 2, "Open questions on assay
priority order") for the open questions this leaves, and
`docs/variant_classification.md` for the current recommended
deduplication approach, which for most categories doesn't consult this
list at all.
"""

ASSAY_PRIORITY_LIST = [
    "BRCA1_Findlay_2018",
    "BRCA2_Hu_2024",
    "VHL_Buckley_2024",
    "JAG1_Gilbert_2024",
    "BARD1_IGVF",
    "PALB2_IGVF",
    "SFPQ_IGVF",
    "RAD51D_IGVF",
    "CTCF_IGVF",
    "BAP1_Waters_2024",
    "DDX3X_Radford_2023",
    "RHO_Wan_2019",
    "RAD51C_Olvera-León_2024",
    "FKRP_Ma_2024",
    "LARGE1_Ma_2024",
    "CARD11_Meitlis_2020_SGE_Ibrutinib_GoF",
    "CARD11_Meitlis_2020_SGE_LoF",
    "ASPA_Grønbæk-Thygesen_2024_abundance",
    "ASPA_Grønbæk-Thygesen_2024_toxicity",
    "BRCA1_Adamovich_2022_Cisplatin_Resistance",
    "BRCA1_Adamovich_2022_HDR",
    "CHEK2_Gebbia_2024",
    "CRX_Shepherdson_2024",
    "F9_Popp_2025_model",
    "G6PD_IGVF",
    "GCK_Gersing_2023_complementation",
    "GCK_Gersing_2024_abundance",
    "KCNE1_Muhammad_2024_trafficking",
    "KCNE1_Muhammad_2024_potassium_flux",
    "KCNE1_Muhammad_2024_trafficking_WT_background_DN",
    "KCNH2_Jiang_2022",
    "KCNH2_Kozek_Glazer_2020",
    "KCNH2_O_Neill_2024_surface_expression",
    "MSH2_Jia_2021",
    "NDUFAF6_Sung_2024",
    "OTC_Lo_2023",
    "PTEN_Matreyek_2018",
    "PTEN_Mighell_2018",
    "SCN5A_Glazer_2020",
    "SCN5A_Ma_2024",
    "SGCB_Li_2023",
    "TP53_Fayer_2021_meta",
    "TP53_Fortuno_2021",
    "TSC2_IGVF",
    "KCNQ4_Zheng_2022_current_homozygous",
    "KCNQ4_Zheng_2022_v12_homozygous",
]
