"""Deduplication priority order for MAVE assays that share a gene.

`ASSAY_PRIORITY_LIST` is used by the OddsPath calibration pipeline
(`notebooks/analysis/OddsPath_classifications.ipynb` and
`notebooks/analysis/Variant_Classification_analysis.ipynb`) to resolve variants that
have functional scores from more than one assay covering the same gene.
Each call site builds `{name: i for i, name in enumerate(ASSAY_PRIORITY_LIST)}`,
maps it onto a `Dataset` column (assays not in the list fall back to
`9999`), then keeps only the lowest-priority-index row per variant via
`sort_values("assay_priority").drop_duplicates(..., keep="first")`.

The list is split into three sections reflecting how each assay's position
was determined:

1. Matches the priority order used in the original preprint submission.
2. Assays that existed at the time of the preprint submission but were not
   part of its explicit priority order. Order is a best guess (grouped by
   gene) and has not been reviewed the way section 1 was.
3. Assays added after the preprint submission. Also a best guess, grouped
   by gene, and unreviewed.

Known open questions: BRCA2_Hu_2024 was and remains prioritized over BRCA2_IGVF
and the Sahu assays; the relative order of the three LDLR_Tabet_2025 assays has
not been decided; and order was never set for genes CBS, HMBS, or PAX6.
"""

ASSAY_PRIORITY_LIST = [
    # --- Section 1: matches the preprint submission pipeline ---
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
    # --- Section 2: existed before the preprint submission but were not ---
    # --- part of its priority order. Best guess, grouped by gene,     ---
    # --- unreviewed.                                                  ---
    "BRCA2_IGVF",
    "BRCA2_Sahu_2025_SGE",
    "BRCA2_Sahu_2023_exon13_SGE",
    "BRCA2_Sahu_2023_exon13_global_score",
    "BRCA2_Sahu_2023_exon13_Cisplatin_Resistance",
    "BRCA2_Sahu_2023_exon13_Olaparib_Resistance",
    "CALM1_CALM2_CALM3_Weile_2017",
    "CBS_Sun_2020_high_B6",
    "CBS_Sun_2020_low_B6",
    "F9_Popp_2025_heavy_chain",
    "F9_Popp_2025_light_chain",
    "F9_Popp_2025_carboxy_F9_specific",
    "F9_Popp_2025_carboxy_gla_motif",
    "F9_Popp_2025_strep_2",
    "HMBS_van_Loggerenberg_2023_combined",
    "HMBS_van_Loggerenberg_2023_erythroid",
    "HMBS_van_Loggerenberg_2023_ubquitous",
    "PAX6_McDonnell_2024_BLX_geneticin",
    "PAX6_McDonnell_2024_BLX_no_geneticin",
    "PAX6_McDonnell_2024_LE9_geneticin",
    "PAX6_McDonnell_2024_LE9_no_geneticin",
    "TARDBP_Bolognesi_Faure_2019",
    "TP53_Kato_2003_AIP1nWT",
    "TP53_Kato_2003_BAXnWT",
    "TP53_Kato_2003_GADD45nWT",
    "TP53_Kato_2003_h1433snWT",
    "TP53_Kato_2003_MDM2nWT",
    "TP53_Kato_2003_NOXAnWT",
    "TP53_Kato_2003_P53R2nWT",
    "TP53_Kato_2003_WAF1nWT",
    "TP53_Giacomelli_2018_combined_score",
    "TP53_Giacomelli_2018_p53WT_Nutlin3",
    "TP53_Giacomelli_2018_p53null_Nutlin3",
    "TP53_Giacomelli_2018_p53null_etoposide",
    "TP53_Boettcher_2019",
    "TPK1_Weile_2017",
    "XRCC2_IGVF",
    # --- Section 3: added after the preprint submission. Best guess,  ---
    # --- grouped by gene, unreviewed.                                 ---
    "BRCA2_Huang_2025_SGE",
    "CHEK2_McCarthy-Leo_2024",
    "LDLR_Tabet_2025_uptake",
    "LDLR_Tabet_2025_abundance",
    "LDLR_Tabet_2025_presence_VLDL",
    "PALB2_Boonen_2026",
    "PALB2_Boonen_2026_SGE",
    "TP53_Funk_2025",
]
