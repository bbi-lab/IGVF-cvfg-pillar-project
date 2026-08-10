import pandas as pd
import numpy as np
from pathlib import Path
from fire import Fire
from functools import reduce
import logging
from io import StringIO
from tqdm import tqdm
import json

logging.basicConfig()
logging.root.setLevel(logging.NOTSET)
logger = logging.getLogger(__name__)


class PillarProjectDataframe:
    def __init__(self, data_path: Path | str):
        self.data_path = Path(data_path)
        self.init_data()

    def init_data(self):
        if not self.data_path.exists():
            raise FileNotFoundError(f"File not found: {self.data_path}")
        self.dataframe = pd.read_csv(self.data_path, sep="\t")

    def __len__(self):
        return len(self.dataframe)

    def get_unique_clinsigs(self):
        sig_sets = self.dataframe.clinvar_sig.apply(
            lambda li: set(_clean_clinsigs(_tolist(li)))
        ).values
        return reduce(lambda x, y: x.union(y), sig_sets)


def _tolist(value, sep="^"):
    try:
        return value.split(sep)
    except AttributeError:
        if pd.isna(value):
            return [
                np.nan,
            ]
        return [
            value,
        ]


def _clean_clinsigs(values):
    return [v.split(";")[0] if isinstance(v, str) else "nan" for v in values]


func_class_map = {#'BRCA1_Findlay_2018': {'LOF':'Abnormal','FUNC':'Normal','INT':'Indeterminate'},
         'BAP1_Waters_2024': {'depleted':'Abnormal','unchanged':'Normal','enriched':'Not specified'},
         'BRCA2_Hu_2024': {'Abnormal': 'Abnormal', 'Normal': 'Normal', 'Intermediate': 'Indeterminate'},
         'CRX_Shepherdson_2024': {'low_activity': 'Abnormal','non-significant': 'Normal',
                                  'high_activity': 'Not specified'},
         'DDX3X_Radford_2023': {'fast depleting': 'Abnormal', 'slow depleting': 'Abnormal',
                                          'unchanged': 'Normal', 'enriched': 'Not specified'},
         'FKRP_Ma_2024': {'damaging_severe': 'Abnormal', 'damaging_mild': 'Abnormal', 
                          'damaging_intermediate': 'Abnormal', 'functional': 'Normal'},
         'JAG1_Gilbert_2024': {'Abnormal': 'Abnormal', 'Likely abnormal': 'Abnormal','Normal': 'Normal'},
         'KCNE1_Muhammad_2024_trafficking_WT_background_DN': {'Loss': 'Abnormal','Possible':'Abnormal','Partial':'Abnormal',
                                               'Normal':'Normal','Gain':'Not specified','PossibleGain':'Not specified'},
         'KCNE1_Muhammad_2024_trafficking': {'Loss': 'Abnormal','Possible':'Abnormal','Partial':'Abnormal',
                                               'Normal':'Normal','Gain':'Not specified','PossibleGain':'Not specified'},
         'KCNE1_Muhammad_2024_potassium_flux': {'Loss': 'Abnormal','Possible':'Abnormal','Partial':'Abnormal',
                                               'Normal':'Normal','Gain':'Not specified','PossibleGain':'Not specified'},
         'LARGE1_Ma_2024': {'damaging':'Abnormal','functional': 'Normal'},
         'NDUFAF6_Sung_2024': {'abnormal': 'Abnormal','normal':'Normal','uncertain': 'Indeterminate'},
         'OTC_Lo_2023': {'Amorphic':'Abnormal','Unimpaired':'Normal', 'Hypomorphic':'Not specified'},
         'RAD51C_Olvera-León_2024': {'fast depleted': 'Abnormal','slow depleted': 'Abnormal',
                                                   'unchanged': 'Normal','enriched':'Not specified'},
         'RHO_Wan_2019':{'low': 'Abnormal','very low': 'Abnormal','high': 'Normal','indeterminate': 'Indeterminate'},
         'SCN5A_Glazer_2020': {'LOF':'Abnormal','possiblyLOF':'Abnormal','possiblyWT':'Normal','WT':'Normal',
                              'GOF': 'Not specified','possiblyGOF':'Not specified'},
         'SCN5A_Ma_2024':{'severe LOF': 'Abnormal', 'moderate LOF': 'Abnormal','normal': 'Normal'},
         'SGCB_Li_2023': {'Non-Functional': 'Abnormal','Functional': 'Normal'},
         'TP53_Fayer_2021_meta': {'Functionally abnormal': 'Abnormal','Functionally normal': 'Normal'},
         'VHL_Buckley_2024': {'LOF1': 'Abnormal','LOF2': 'Abnormal','Neutral': 'Normal', 'Intermediate': 'Indeterminate'},
         'BARD1_IGVF': {'functionally_abnormal': 'Abnormal', 'functionally_normal': 'Normal', 'indeterminate': 'Indeterminate'},
         'PALB2_IGVF': {'functionally_abnormal': 'Abnormal', 'functionally_normal': 'Normal', 'indeterminate': 'Indeterminate'},
         'CTCF_IGVF': {'functionally_abnormal': 'Abnormal', 'functionally_normal': 'Normal', 'indeterminate': 'Indeterminate'},
         'RAD51D_IGVF': {'functionally_abnormal': 'Abnormal', 'functionally_normal': 'Normal', 'indeterminate': 'Indeterminate'},
         'SFPQ_IGVF': {'functionally_abnormal': 'Abnormal', 'functionally_normal': 'Normal', 'indeterminate': 'Indeterminate'},
         'PTEN_Matreyek_2018': {'functionally_abnormal': 'Abnormal', 'functionally_normal': 'Normal'},
         'F9_Popp_2025_model': {'WT-like': 'Normal','Loss of function': 'Abnormal'},
         'G6PD_IGVF': {'functionally_abnormal': 'Abnormal', 'functionally_normal': 'Normal'},
         'TSC2_IGVF': {'functionally_abnormal': 'Abnormal', 'functionally_normal': 'Normal'},
         'CARD11_Meitlis_2020_SGE_LoF': {'functional': 'Normal', 'not definitive': 'Indeterminate',
                                                 'likely functional': 'Normal','likely nonfunctional': 'Abnormal', 
                                                 'nonfunctional': 'Abnormal'}}

INTERVAL_CLASS_NORMALIZATION = {
    "normal": "Normal",
    "abnormal": "Abnormal",
    "not specified": "Indeterminate",
}

def classify_by_intervals(score, intervals):
    """
    intervals: list of dicts with keys:
        - 'range': tuple (start, end, start_inclusive, end_inclusive)
        - 'class': raw MaveDB class string
    """
    for interval in intervals:
        lo, hi, lo_inc, hi_inc = interval["range"]

        in_range = (
            (score > lo or (lo_inc and score == lo)) and
            (score < hi or (hi_inc and score == hi))
        )

        if in_range:
            raw = interval["class"]
            if raw is None or str(raw).lower() == "nan":
                return "Indeterminate"

            return INTERVAL_CLASS_NORMALIZATION.get(
                raw.strip().lower(),
                "Indeterminate"
            )

    return "Indeterminate"

import math
import re

def parse_interval_range(range_str):
    """
    Parses strings like:
      '[-0.748, Inf)'
      '(-Inf,-1.328]'
    Returns: (lo, hi, lo_inc, hi_inc)
    """
    if not isinstance(range_str, str):
        return None

    m = re.match(r'([\[\(])\s*([^,]+)\s*,\s*([^\]\)]+)\s*([\]\)])', range_str)
    if not m:
        return None

    lo_br, lo, hi, hi_br = m.groups()

    def to_float(x):
        x = x.strip()
        if x.lower() in {"-inf", "-infinity"}:
            return -math.inf
        if x.lower() in {"inf", "infinity"}:
            return math.inf
        return float(x)

    return (
        to_float(lo),
        to_float(hi),
        lo_br == "[",
        hi_br == "]",
    )

def extract_intervals_from_row(row, max_intervals=6):
    intervals = []

    for i in range(1, max_intervals + 1):
        range_key = f"Interval {i} range"
        class_key = f"Interval {i} MaveDB class"

        range_str = row.get(range_key)
        class_val = row.get(class_key)

        parsed = parse_interval_range(range_str)
        if parsed is None:
            continue

        intervals.append({
            "range": parsed,
            "class": class_val
        })

    return intervals



def standardize_auth_label(
    dataset_name,
    func_class,
    score=None,
    intervals=None,
):
    """
    Returns standardized label: Normal / Abnormal / Indeterminate

    Priority:
      1) Dataset-specific func_class_map
      2) Interval-based classification using score
      3) Indeterminate
    """

    if dataset_name in func_class_map:
        dataset_map = func_class_map[dataset_name]
        if func_class in dataset_map:
            return dataset_map[func_class]

    if isinstance(func_class, str):
        norm = func_class.strip().lower()
        if norm in INTERVAL_CLASS_NORMALIZATION:
            return INTERVAL_CLASS_NORMALIZATION[norm]

    if score is not None and intervals is not None:
        return classify_by_intervals(score, intervals)

    return "Indeterminate"



class BasicScoreset:
    def __init__(self, scores: np.ndarray, sample_assignments: np.ndarray,**kwargs):
        self.scores = scores
        self._sample_assignments = sample_assignments
        self.validate_inputs()
        self.validate_sample_assignments()
        self.sample_counts = self._sample_assignments.sum(axis=0)
        self.scoreset_name = kwargs.get("scoreset_name", "BasicScoreset")
        self.sample_names = kwargs.get("sample_names",["Pathogenic/Likely Pathogenic","Benign/Likely Benign","gnomAD","Synonymous"])[:self.n_samples]

    def validate_inputs(self):
        n_observations = self.scores.shape[0]
        if self._sample_assignments.shape[0] != n_observations:
            raise ValueError(
                f"Number of observations in scores {n_observations:,d} does not match number of rows in sample_assignments {self._sample_assignments.shape[0]:,d}"
            )

    def validate_sample_assignments(self):
        ndim = self._sample_assignments.ndim
        if ndim == 1:
            print(
                "Assuming sample_assignments is a list of sample identifiers, converting to 2D array."
            )
            sample_ids = np.array(self._sample_assignments)
            unique_samples = list(set(sample_ids))
            self._sample_assignments = np.zeros(
                (len(sample_ids), len(unique_samples)), dtype=bool
            )
            for sampleNum, sample_id in enumerate(unique_samples):
                self._sample_assignments[:, sampleNum] = sample_ids == sample_id
        elif ndim != 2:
            raise ValueError(
                f"sample_assignments must be a 1D list of sample ids or 2D array of one-hot vectors, got {ndim} dimensions"
            )

    @property
    def sample_assignments(self):
        return self._sample_assignments[:, self.sample_counts > 0]

    @property
    def n_samples(self):
        return self._sample_assignments.shape[1]

    @property
    def samples(self):
        for sample_index in range(self.n_samples):
            if self.sample_counts[sample_index] > 0:
                yield self.scores[
                    self._sample_assignments[:, sample_index]
                ], self.sample_names[sample_index]


    # @property
    # def samples(self):
    #     """
    #     Iterate over the samples in the scoreset, yielding the sample scores and sample name.
    #     """
    #     for sample_index in range(self.sample_assignments.shape[1]):
    #         sample_scores = self.scores[self.sample_assignments[:, sample_index]]
    #         if len(sample_scores) > 0:
    #             yield sample_scores, f"Sample {sample_index + 1}"

    @classmethod
    def from_csv(cls, csv_path: Path | str, **kwargs):
        """
        Create a BasicScoreset from a CSV file.

        Required columns:
         - scores: The scores for each observation
         - sample_assignments: The ID of the sample to which each observation belongs
            Example:
            scores,sample_assignments
            0.5,1
            0.7,2
            0.3,1

        Parameters
        ----------
        csv_path : Path|str
            The path to the CSV file to create the scoreset from
        """
        csv_path = Path(csv_path)
        if not csv_path.exists():
            raise FileNotFoundError(f"CSV file not found: {csv_path}")
        df = pd.read_csv(csv_path)
        if "scores" not in df.columns or "sample_assignments" not in df.columns:
            raise ValueError(
                "CSV must contain 'scores' and 'sample_assignments' columns"
            )
        scores = np.array(df["scores"].values)
        sample_assignments = np.array(df["sample_assignments"].values)
        return cls(scores, sample_assignments, **kwargs)


class Scoreset:
    def __init__(self, dataframe: pd.DataFrame, **kwargs):
        """
        Required Arg:
        - dataframe : pd.DataFrame : Pillar Project dataframe

        Optional Arg:
        - min_clinvar_star : int (default 1) : minimum review status (star count) to use Clinical Significance annotations
        - clinvar_release : str in {'2025','2018'} (default '2025') : Clinvar release to use
        """
        self._init_kwargs = dict(kwargs)
        self._init_dataframe(dataframe, **kwargs)

    def to_json(self, output_path: Path | str):
        """
        Save the scoreset to a JSON file.

        Parameters
        ----------
        output_path : Path|str
            The path to save the JSON file to

        Returns
        -------
        None
        """
        output_path = Path(output_path)
        if not output_path.parent.exists():
            output_path.parent.mkdir(parents=True, exist_ok=True)
        self.dataframe.to_json(output_path, orient="records", lines=True)

    @classmethod
    def from_dataframe(cls, dataframe: pd.DataFrame, **kwargs):
        """
        Create a Scoreset from a pandas DataFrame.

        Parameters
        ----------
        dataframe : pd.DataFrame
            The dataframe to create the scoreset from

        Returns
        -------
        Scoreset
            A Scoreset object initialized with the given dataframe
        """
        return cls(dataframe, **kwargs)

    @classmethod
    def from_json(cls, json_path: Path | str, **kwargs):
        """
        Create a Scoreset from a JSON or JSONL file.
    
        Automatically detects:
          - Standard JSON object or list
          - Newline-delimited JSON (JSONL)
    
        Parameters
        ----------
        json_path : Path|str
            Path to JSON/JSONL input file
    
        Returns
        -------
        Scoreset
            A Scoreset initialized with parsed data
        """
        json_path = Path(json_path)
        if not json_path.exists():
            raise FileNotFoundError(f"JSON file not found: {json_path}")
    
        with open(json_path, "r") as f:
            text = f.read().strip()
    
        # Try standard JSON first
        try:
            data = json.loads(text)
        except json.JSONDecodeError as e:
            # Only fallback if this specific error indicates multiple objects
            if "Extra data" not in str(e):
                raise  # If it's a different JSON error, report it normally
    
            # Fallback for JSONL: parse line-by-line
            try:
                lines = [json.loads(line) for line in text.splitlines() if line.strip()]
                data = lines  # list of JSON objects
            except Exception:
                raise ValueError(
                    f"File is neither valid JSON nor JSONL: {json_path}\n"
                    f"Original error: {e}"
                )
    
        # Convert result into DataFrame
        dataframe = pd.DataFrame.from_records(data)
        return cls(dataframe, **kwargs)


    def _init_dataframe(self, dataframe: pd.DataFrame, **kwargs):
        """
        Initialize the scoreset from the dataframe

        Parameters
        ----------
        dataframe : pd.DataFrame
            The dataframe to initialize the scoreset from

        Optional Arg:
            - min_clinvar_star : int (default 1) : minimum review status (star count) to use Clinical Significance annotations
            - clinvar_release : str in {'2025','2018'} (default '2025') : Clinvar release to use
        Returns
        -------
        None
        """
        if not isinstance(dataframe, pd.DataFrame):
            raise TypeError("dataframe must be a pandas DataFrame")
        if len(dataframe.Dataset.unique()) != 1:
            raise ValueError("dataframe must contain only one dataset")
        if not len(dataframe):
            raise ValueError("dataframe must contain at least one row")
        # drop rows with NaN in auth_reported_score
        dataframe = dataframe.assign(
            auth_reported_score=pd.to_numeric(
                dataframe.auth_reported_score, errors="coerce"
            )
        )
        dataframe = dataframe.dropna(subset=["auth_reported_score"])
        dataframe = Scoreset.remove_outliers(dataframe, **kwargs)
        if not len(dataframe):
            raise ValueError(
                "dataframe must contain at least one row with a non-NaN auth_reported_score"
            )
        self.dataframe = dataframe
        self.filter_nonsense = kwargs.get("filter_nonsense",False)
        self.filter_invalid()
        self.splicing_filter(**kwargs)
        min_clinvar_star = kwargs.get("min_clinvar_star",1)
        clinvar_release = kwargs.get("clinvar_release",'2025')
        self.variants = [Variant(row,min_clinvar_star,clinvar_release) for _, row in self.dataframe.iterrows()]
        self.synonymous_exclusive = kwargs.get("synonymous_exclusive",True)
        self.score_avg = kwargs.get("score_avg",False)
        self._init_matrices(**kwargs)

    def filter_invalid(self):
        if "Flag" in self.dataframe.columns:
            self.dataframe = self.dataframe[self.dataframe.Flag != "*"]

        if self.filter_nonsense:
            # print('filtering nonsense variants within 50 aa of termini...')
            mask = (
                # Keep all non-stop_gained variants
                (self.dataframe['simplified_consequence'] != 'stop_gained') |
                # OR keep stop_gained only if position is between 51-349
                ((self.dataframe['simplified_consequence'] == 'stop_gained') &
                 (self.dataframe['aa_pos'] >= 51) &
                 (self.dataframe['aa_pos'] <= 349))
            )
            self.dataframe = self.dataframe[mask]

    def splicing_filter(self, **kwargs):
        self.detects_splice = (
            self.dataframe.loc[:, "splice_measure"].unique()[0] == "Yes"  # type: ignore
        )
        if not self.detects_splice:
            self.dataframe = self.dataframe[
                self.dataframe.simplified_consequence.str.lower() != "splice region"
            ]
            self.dataframe = self.dataframe[
                self.dataframe.simplified_consequence.str.lower() != "splice_site_variant"
            ]
            self.dataframe = self.dataframe[
                (self.dataframe["spliceAI_DS_AG"] < 0.2) &
                (self.dataframe["spliceAI_DS_AL"] < 0.2) &
                (self.dataframe["spliceAI_DS_DG"] < 0.2) &
                (self.dataframe["spliceAI_DS_DL"] < 0.2)
            ]

    @staticmethod
    def remove_outliers(dataframe, **kwargs):
        """
        Optionally clip the dataframe to remove observations outside a specified percentile range

        Parameters
        ----------
        dataframe : pd.DataFrame
            The dataframe to remove outliers from

        Optional Parameters
        -------------------
        - quantile_min : float (default: 0.0)
        - quantile_max : float (default: 1.0)

        Returns
        -------
        pd.DataFrame
            The dataframe with outliers removed (1.5 IQR Rule)
        """
        quantile_min = kwargs.get("quantile_min", 0.0)
        quantile_max = kwargs.get("quantile_max", 1.0)
        lowerbound = dataframe.auth_reported_score.quantile(quantile_min)
        upperbound = dataframe.auth_reported_score.quantile(quantile_max)
        scores = dataframe.auth_reported_score
        include = (scores >= lowerbound) & (scores <= upperbound)
        return dataframe[include]

    def __len__(self):
        return len(self.variants)

    def _init_matrices(self, **kwargs):
        """
        Optional Arguments:
        - population_type : str : one of {'all_variants',
                                  'all_nsSNV',
                                  'all_missense_nsSNV',
                                  'gnomAD',
                                  'gnomAD_nsSNV',
                                  'gnomAD_missense_nsSNV'}
        
        """
        self.has_synomyous = any([variant.is_synonymous for variant in self.variants])
        # if self.has_synomyous:
        self.NSamples = 4
        self.sample_names = [
            "Pathogenic/Likely Pathogenic",
            "Benign/Likely Benign",
            "population",
            "Synonymous",
        ]
        # else:
        #     self.NSamples = 3
        #     self.sample_names = [
        #         "Pathogenic/Likely Pathogenic",
        #         "Benign/Likely Benign",
        #         "population",
        #     ]
        variants_by_id = self.get_variants_by_id()
        self.n_variants = len(variants_by_id)
        self._sample_assignments = np.zeros(
            (self.n_variants, self.NSamples), dtype=bool
        )
        self._scores = np.zeros(self.n_variants)
        self._auth_labels = np.empty(self.n_variants, dtype='U50')
        self._snv_scores = []
        self._vus_scores = []
        self._all_scores = []
        self._aa_subs = np.empty(self.n_variants, dtype='U50')
        self._ids = []
        self.parse_population_type(**kwargs)
        for idx, (_id, variants) in enumerate(variants_by_id.items()):

                
            self._ids.append(_id)

            if self.score_avg:
                variants[0].auth_reported_score = np.mean([v.auth_reported_score for v in variants])
            
            self._scores[idx] = variants[0].auth_reported_score
            self._auth_labels[idx] = variants[0].auth_label
            if variants[0].is_snv:
                self._snv_scores.append(variants[0].auth_reported_score)
            self._all_scores.append(variants[0].auth_reported_score)
            
            if any([variant.is_vus for variant in variants]):
                self._vus_scores.append(variants[0].auth_reported_score)
            
            try:
                self._aa_subs[idx] = variants[0].aa_ref + str(int(variants[0].aa_pos)) + variants[0].aa_alt
            except (ValueError, TypeError):
                self._aa_subs[idx] = None
            if any([variant.is_synonymous for variant in variants]):
                self._sample_assignments[idx, 3] = True
                if self.synonymous_exclusive:
                    continue
            if any([self.is_population_member(variant) for variant in variants]):
                self._sample_assignments[idx, 2] = True
            if any([variant.is_pathogenic for variant in variants]):
                self._sample_assignments[idx, 0] = True
            if any([variant.is_benign for variant in variants]):
                self._sample_assignments[idx, 1] = True
        
        keep_mask = self._sample_assignments.any(axis=1)
        self._keep_mask = keep_mask
        self._scores = self.scores[keep_mask]
        self._sample_assignments = self._sample_assignments[keep_mask]
        self._auth_labels = self._auth_labels[keep_mask]
        self.n_variants = len(self._scores)
        self.sample_counts = self._sample_assignments.sum(axis=0)
        self._snv_scores = np.array(self._snv_scores)
        self._vus_scores = np.array(self._vus_scores)
        self._all_scores = np.array(self._all_scores)

    def parse_population_type(self,**kwargs):
        population_type = kwargs.get("population_type",'gnomAD')
        valid_population_types = {'all_variants',
                                  'all_nsSNV',
                                  'all_missense_nsSNV',
                                  'gnomAD',
                                  'gnomAD_nsSNV',
                                  'gnomAD_missense_nsSNV'}
        if population_type not in valid_population_types:
            raise ValueError(f"Invalid population type {population_type}; must be in {valid_population_types}")
        self.population_type = population_type

    def is_population_member(self,variant)->bool:
        if self.population_type == "all_variants":
            return True
        if self.population_type == "all_nsSNV":
            return variant.is_nsSNV
        if self.population_type == "all_missense_nsSNV":
            return variant.simplified_consequence == "missense_variant" and variant.is_nsSNV
        if self.population_type == "gnomAD":
            return variant.is_gnomAD
        if self.population_type == "gnomAD_nsSNV":
            return variant.is_gnomAD and variant.is_nsSNV
        if self.population_type == "gnomAD_missense_nsSNV":
            return variant.is_gnomAD and \
                variant.is_nsSNV and \
                variant.simplified_consequence == "missense_variant"
        return False
        

    def get_variants_by_id(self):
        """
        Iterate over all unique Variant.ID values, returning the variants with that given ID.

        Returns
        -------
        dict
            A dictionary where keys are unique Variant.ID values and values are lists of Variant objects with that ID
        """
        variants_by_id = {}
        for variant_id in set(variant.ID for variant in self.variants):
            variants_by_id[variant_id] = [
                variant for variant in self.variants if variant.ID == variant_id
            ]
        return variants_by_id

    @property
    def sample_assignments(self):
        return self._sample_assignments[:, self.sample_counts > 0]

    @property
    def n_samples(self):
        return self.sample_assignments.shape[1]

    @property
    def samples(self):
        for sample_index in range(self.NSamples):
            if self.sample_counts[sample_index] > 0:
                yield self.scores[
                    self._sample_assignments[:, sample_index]
                ], self.sample_names[sample_index]

    @property
    def scores(self):
        return self._scores

    @property
    def snv_scores(self):
        return self._snv_scores
        
    @property
    def vus_scores(self):
        return self._vus_scores
    
    @property
    def all_scores(self):
        return self._all_scores
        
    @property
    def aa_subs(self):
        return self._aa_subs

    @property
    def auth_labels(self):
        return self._auth_labels

    @auth_labels.setter
    def auth_labels(self, value):
        self._auth_labels = value

    @property
    def scoreset_name(self):
        return self.dataframe.Dataset.values[0]

    @property
    def keep_mask(self):
        return self._keep_mask

    def __repr__(self):
        out = f"{self.scoreset_name}: {len(self)} total variants\n"
        for sample_scores, sample_name in self.samples:
            out += f"\t{sample_name}: {len(sample_scores)} variants\n"

        return out

    def to_serializable(self):
        return {
            "__type__": "Scoreset",
            "dataframe": self.dataframe.to_dict(orient="split"),
            "kwargs": self._init_kwargs,
        }

    
    @classmethod
    def from_serializable(cls, payload: dict):
        df = pd.DataFrame(**payload["dataframe"])
        return cls(df, **payload["kwargs"])



class Variant:
    def __init__(self, variant_info: pd.Series, min_clinvar_star: int, clinvar_release: str):
        self.min_clinvar_star = min_clinvar_star
        self.clinvar_release = clinvar_release
        self.row = variant_info
        self._init_variant_info(variant_info)

    def _init_variant_info(self, variant_info: pd.Series):
        self.ID = None
        self.simplified_consequence = None
        self.clinvar_star = None
        self.gnomad_MAF = None
        self.auth_reported_score = None
        self.transcript_ref = ""
        self.transcript_alt = ""
        self.aa_ref = np.nan
        self.aa_alt = ""
        self.hgvs_p = ""
        for k, v in variant_info.items():
            setattr(self, str(k), v)
        self.clinvar_sig = getattr(variant_info,f"clinvar_sig_{self.clinvar_release}")
        self.clinvar_star = getattr(variant_info,f"clinvar_star_{self.clinvar_release}")
        self.parse_gnomAD_MAF()
        self.parse_clinvar_sig()
        self.parse_consequences()
        self.assign_nsSNV()
        self.parse_auth_class()

    def assign_nsSNV(self):
        # is the variant a nonsynonymous nsSNV
        self.is_nsSNV = (isinstance(self.transcript_ref,str) and len(self.transcript_ref) == 1) and \
                        (isinstance(self.transcript_alt,str) and len(self.transcript_alt) == 1) and \
                        (isinstance(self.hgvs_p,float) | \
                         (self.aa_ref != self.aa_alt))

    def parse_consequences(self):
        self.is_synonymous = (self.simplified_consequence == "Synonymous") or (
            self.simplified_consequence == "synonymous_variant"
        )
        self.is_missense = self.simplified_consequence == 'missense_variant'
        self.is_snv = len(str(self.ref_allele)) == 1 and len(str(self.alt_allele)) == 1 and str(self.ref_allele) != str(self.alt_allele)

    def parse_clinvar_sig(self):
        self.is_conflicting = (
            self.clinvar_sig == "Conflicting classifications of pathogenicity"
        )
        self.eval_quality()
        
        self.is_benign = self.sufficient_quality and self.clinvar_sig in {
            "Benign",
            "Likely benign",
            "Benign/Likely benign",
        }
        self.is_pathogenic = self.sufficient_quality and self.clinvar_sig in {
            "Pathogenic",
            "Likely pathogenic",
            "Pathogenic/Likely pathogenic",
        }
        self.is_vus = self.sufficient_quality_2025 and self.clinvar_sig_2025 in { # VUS ALWAYS 2025 SINCE NOT CONTROL
            "Uncertain significance",
        }

    def standardize_class(self, c):
        # c: Class
        if c is None:
            return c
        elif pd.isna(c) or (isinstance(c, str) and c.strip().lower() in {"nan", "not specified"}):
            return "Indeterminate"
        else:
            return c.capitalize()
        

    def parse_auth_class(self):
        self.auth_label = self.standardize_class(getattr(self, "StandardizedClass", None))
        if self.auth_label is not None:
            return
    
        func_class = getattr(self, "auth_reported_func_class", None)
    
        score = getattr(self, "score", None)
    
        intervals = extract_intervals_from_row(self.row)
        # print('intervals',intervals)
    
        if func_class is not None and self.Dataset in func_class_map:
            self.auth_label = standardize_auth_label(
                dataset_name=self.Dataset,
                func_class=func_class,
                score=score,
                intervals=intervals if intervals else None
            )
            # print("class",self.auth_label)
        else:
            if score is not None and intervals:
                self.auth_label = classify_by_intervals(score, intervals)
                # print("interval",self.auth_label)
            else:
                self.auth_label = "Indeterminate"
                # print("ind",self.auth_label)


    def eval_quality(self):
        one_star_status = {
            'criteria provided, single submitter',
            'criteria provided, conflicting interpretations',
            'criteria provided, conflicting classifications',

        }
        two_star_statuses = {
            'criteria provided, multiple submitters, no conflicts',

        }
        three_star_statuses = {
            'reviewed by expert panel'
        }
        zero_star_statuses = {
            np.nan,
            'no assertion criteria provided',
            'no assertion provided',
            'no interpretation for the single variant',
            'no classification provided',
            '-'
        }
        if self.min_clinvar_star == 0:
            self.sufficient_quality = True
            self.sufficient_quality_2025 = True
        elif self.min_clinvar_star == 1:
            self.sufficient_quality = self.clinvar_star not in zero_star_statuses
            self.sufficient_quality_2025 = self.clinvar_star_2025 not in zero_star_statuses
        elif self.min_clinvar_star == 2:
            self.sufficient_quality = self.clinvar_star not in zero_star_statuses.union(one_star_status)
            self.sufficient_quality_2025 = self.clinvar_star_2025 not in zero_star_statuses.union(one_star_status)
        elif self.min_clinvar_star == 3:
            self.sufficient_quality = self.clinvar_star not in zero_star_statuses.union(one_star_status).union(two_star_statuses)
            self.sufficient_quality_2025 = self.clinvar_star_2025 not in zero_star_statuses.union(one_star_status).union(two_star_statuses)
        elif self.min_clinvar_star == 4:
            self.sufficient_quality = self.clinvar_star not in zero_star_statuses.union(one_star_status).union(two_star_statuses).union(three_star_statuses)
            self.sufficient_quality_2025 = self.clinvar_star_2025 not in zero_star_statuses.union(one_star_status).union(two_star_statuses).union(three_star_statuses)
        else:
            raise ValueError(f"Invalid min_clinvar_star value {self.min_clinvar_star}")

    def parse_gnomAD_MAF(self):
        """
        It is possible that the MAF is a list of values separated by a semicolon. If so, parse the list and obtain the maximum value.
        """
        self.is_gnomAD = not pd.isna(self.gnomad_MAF)

    @property
    def score(self):
        return self.auth_reported_score

    @staticmethod
    def is_nan(value):
        return pd.isna(value) or value == "nan"


def summarize_datasets(dataframe_path, **kwargs):
    """
    Summarize the datasets in the dataframe at dataframe_path.

    Parameters
    ----------
    dataframe_path : str
        The path to the dataframe containing the dataset

    Keyword Arguments
    -----------------
    - output_file : str|Path
        The path to save the summary to

    Returns
    -------
    None
    """
    output_file = kwargs.get("output_file", None)
    if output_file is not None:
        output_file = Path(output_file)
        # output_file.mkdir(parents=True, exist_ok=True)
        f = open(str(output_file), "w")
    else:
        f = StringIO()
    df = PillarProjectDataframe(dataframe_path)
    for dataset_name, ds_df in df.dataframe.groupby("Dataset"):
        scoreset = Scoreset(
            ds_df,
            missense_only=kwargs.get("missense_only", False),
            synonymous_exclusive=kwargs.get("synonymous_exclusive", True),
            score_avg=kwargs.get("score_avg", False),
            filter_nonsense=kwargs.get("filter_nonsense", False)
        )
        f.write(f"{dataset_name}\n")
        f.write(str(scoreset))
        f.write("\n")
    if isinstance(f, StringIO):
        print(f.getvalue())
    else:
        f.close()


def csv_to_vcf(input_filepath, output_filepath):
    """
    Convert a CSV file to a gzipped VCF file.

    Parameters
    ----------
    input_filepath : str|Path
        The path to the input CSV file.
    output_filepath : str|Path
        The path to the output gzipped VCF file.

    Returns
    -------
    None
    """
    input_filepath = Path(input_filepath)
    output_filepath = Path(output_filepath)

    # Read the CSV file
    df = pd.read_csv(input_filepath)

    # Filter rows with non-null hg38_start
    df = df[df["hg38_start"].notnull()]

    # Open the output file for writing
    with open(output_filepath, "w") as vcf_file:
        # Write VCF header
        vcf_file.write("##fileformat=VCFv4.2\n")
        vcf_file.write("##source=tsv_to_vcf\n")
        vcf_file.write(
            """##contig=<ID=1,length=248956422,assembly=GRCh38>
##contig=<ID=2,length=242193529,assembly=GRCh38>
##contig=<ID=3,length=198295559,assembly=GRCh38>
##contig=<ID=4,length=190214555,assembly=GRCh38>
##contig=<ID=5,length=181538259,assembly=GRCh38>
##contig=<ID=6,length=170805979,assembly=GRCh38>
##contig=<ID=7,length=159345973,assembly=GRCh38>
##contig=<ID=8,length=145138636,assembly=GRCh38>
##contig=<ID=9,length=138394717,assembly=GRCh38>
##contig=<ID=10,length=133797422,assembly=GRCh38>
##contig=<ID=11,length=135086622,assembly=GRCh38>
##contig=<ID=12,length=133275309,assembly=GRCh38>
##contig=<ID=13,length=114364328,assembly=GRCh38>
##contig=<ID=14,length=107043718,assembly=GRCh38>
##contig=<ID=15,length=101991189,assembly=GRCh38>
##contig=<ID=16,length=90338345,assembly=GRCh38>
##contig=<ID=17,length=83257441,assembly=GRCh38>
##contig=<ID=18,length=80373285,assembly=GRCh38>
##contig=<ID=19,length=58617616,assembly=GRCh38>
##contig=<ID=20,length=64444167,assembly=GRCh38>
##contig=<ID=21,length=46709983,assembly=GRCh38>
##contig=<ID=22,length=50818468,assembly=GRCh38>
##contig=<ID=X,length=156040895,assembly=GRCh38>
##contig=<ID=Y,length=57227415,assembly=GRCh38>
##contig=<ID=M,length=16569,assembly=GRCh38>
"""
        )
        vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        df = df.sort_values(by=["Chrom", "hg38_start"])  # type: ignore
        # Write VCF rows
        for _, row in tqdm(df.iterrows(), total=len(df)):
            vcf_file.write(
                f"{row['Chrom']}\t{int(row.hg38_start)}\t{row['ID']}\t{row['ref_allele']}\t{row['alt_allele']}\t.\t.\t.\n"
            )


if __name__ == "__main__":
    # summarize_datasets("../../outputs/integrated_variant_effect_dataset.tsv.gz",missense_only=False, synonymous_exclusive=False,output_file="dataset_summary_all_synonymousNonExclusive.txt")
    Fire()
    # summarize_datasets("/data/dzeiberg/pillar_project/dataframe/pillar_data_condensed_gold_standard_02_05_25.csv",missense_only=False, synonymous_exclusive=False,output_file="dataset_summary_all_synonymousNonExclusive.txt")