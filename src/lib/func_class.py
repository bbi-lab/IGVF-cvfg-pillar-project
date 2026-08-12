"""Historical per-dataset author-reported functional class label mapping.

`FUNC_CLASS_LABEL_MAP` was used by `annotate_func_class` in
`notebooks/analysis/Variant_Classification_analysis.ipynb` to translate each dataset's
author-reported functional classification labels (e.g. "depleted",
"loss-of-function", "Amorphic") into this project's standardized functional
class ontology (Normal, Abnormal, Indeterminate).

As of 2026-08, functional class categories are instead taken directly from
MaveDB via `auth_reported_func_class_category`, so classification stays
consistent with MaveDB's own categorization rather than a manually maintained
mapping. This dictionary is no longer used by any pipeline code -- it's kept
here for documentation/historical reference only.
"""

FUNC_CLASS_LABEL_MAP = {
    "BAP1_Waters_2024": {"depleted": "Abnormal", "unchanged": "Normal", "enriched": "Not specified"},
    "BRCA2_Hu_2024": {
        "Functionally abnormal": "Abnormal",
        "Functionally normal": "Normal",
        "Intermediate": "Indeterminate",
    },
    "CRX_Shepherdson_2024": {
        "low_activity": "Abnormal",
        "non-significant": "Normal",
        "high_activity": "Not specified",
    },
    "DDX3X_Radford_2023": {
        "fast depleting": "Abnormal",
        "slow depleting": "Abnormal",
        "unchanged": "Normal",
        "enriched": "Not specified",
    },
    "FKRP_Ma_2024": {
        "damaging severe": "Abnormal",
        "damaging mild": "Abnormal",
        "damaging intermediate": "Abnormal",
        "functional": "Normal",
    },
    "JAG1_Gilbert_2024": {"abnormal": "Abnormal", "likely abnormal": "Abnormal", "normal": "Normal"},
    "KCNE1_Muhammad_2024_trafficking_WT_background_DN": {
        "Loss": "Abnormal",
        "Possible": "Abnormal",
        "Partial": "Abnormal",
        "Normal": "Normal",
        "Gain": "Not specified",
        "PossibleGain": "Not specified",
    },
    "KCNE1_Muhammad_2024_trafficking": {
        "Loss": "Abnormal",
        "Possible": "Abnormal",
        "Partial": "Abnormal",
        "Normal": "Normal",
        "Gain": "Not specified",
        "PossibleGain": "Not specified",
    },
    "KCNE1_Muhammad_2024_potassium_flux": {
        "loss-of-function": "Abnormal",
        "partial loss-of-function": "Abnormal",
        "normal function": "Normal",
        "gain-of-function": "Not specified",
        "PossibleGain": "Not specified",
        # No longer present: 'Possible':'Abnormal', 'PossibleGain':'Not specified'
    },
    "LARGE1_Ma_2024": {"damaging": "Abnormal", "functional": "Normal"},
    "NDUFAF6_Sung_2024": {"abnormal": "Abnormal", "normal": "Normal", "uncertain": "Indeterminate"},
    "OTC_Lo_2023": {
        "Amorphic": "Abnormal",
        "Functionally unimpaired": "Normal",
        "Hypomorphic": "Not specified",
    },
    "RAD51C_Olvera-León_2024": {
        "Fast depleted": "Abnormal",
        "Slow depleted": "Abnormal",
        "Unchanged": "Normal",
        "Enriched": "Not specified",
    },
    "RHO_Wan_2019": {
        "low": "Abnormal",
        "very low": "Abnormal",
        "high": "Normal",
        "indeterminate": "Indeterminate",
    },
    "SCN5A_Glazer_2020": {
        "LOF": "Abnormal",
        "possiblyLOF": "Abnormal",
        "possiblyWT": "Normal",
        "WT": "Normal",
        "GOF": "Not specified",
        "possiblyGOF": "Not specified",
    },
    "SCN5A_Ma_2024": {"severe LOF": "Abnormal", "moderate LOF": "Abnormal", "normal": "Normal"},
    "SGCB_Li_2023": {"Non-functional": "Abnormal", "Functional": "Normal"},
    "TP53_Fayer_2021_meta": {"Functionally abnormal": "Abnormal", "Functionally normal": "Normal"},
    "VHL_Buckley_2024": {
        "LOF1": "Abnormal",
        "LOF2": "Abnormal",
        "Neutral": "Normal",
        "Intermediate": "Indeterminate",
    },
    "BARD1_IGVF": {
        "Functionally Abnormal": "Abnormal",
        "Functionally Normal": "Normal",
        "Indeterminate": "Indeterminate",
    },
    "PALB2_IGVF": {
        "Functionally Abnormal": "Abnormal",
        "Functionally Normal": "Normal",
        "Indeterminate": "Indeterminate",
    },
    "CTCF_IGVF": {
        "Functionally Abnormal": "Abnormal",
        "Functionally Normal": "Normal",
        "Indeterminate": "Indeterminate",
    },
    "RAD51D_IGVF": {
        "Functionally Abnormal": "Abnormal",
        "Functionally Normal": "Normal",
        "Indeterminate": "Indeterminate",
    },
    "SFPQ_IGVF": {
        "Functionally Abnormal": "Abnormal",
        "Functionally Normal": "Normal",
        "Indeterminate": "Indeterminate",
    },
    "XRCC2_IGVF": {
        "Functionally Abnormal": "Abnormal",
        "Functionally Normal": "Normal",
        "Indeterminate": "Indeterminate",
    },
    "PTEN_Matreyek_2018": {
        "low": "Abnormal",
        "possibly_low": "Abnormal",
        "possibly_wt-like": "Normal",
        "wt-like": "Normal",
        # Formerly just {'functionally_abnormal': 'Abnormal', 'functionally_normal': 'Normal'},
    },
    "F9_Popp_2025_model": {"WT-like": "Normal", "Loss of function": "Abnormal"},
    "G6PD_IGVF": {"Functionally Abnormal": "Abnormal", "Functionally Normal": "Normal"},
    "TSC2_IGVF": {"Functionally Abnormal": "Abnormal", "Functionally Normal": "Normal"},
    "CARD11_Meitlis_2020_SGE_LoF": {
        "functional": "Normal",
        "not definitive": "Indeterminate",
        "likely functional": "Normal",
        "likely nonfunctional": "Abnormal",
        "nonfunctional": "Abnormal",
    },
    "CARD11_Meitlis_2020_SGE_Ibrutinib_GoF": {
        "likely not gain of function": "Normal",
        "not definitive": "Indeterminate",
        "not gain of function": "Normal",
        "likely gain of function": "Abnormal",
        "gain of function": "Abnormal",
    },
}
