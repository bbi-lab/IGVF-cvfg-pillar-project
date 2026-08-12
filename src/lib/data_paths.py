"""Shared input/output directory layout for the OddsPath/ExCALIBR pipeline.

`get_data_paths` is used by `notebooks/analysis/OddsPath_calculations.ipynb`,
`notebooks/analysis/Variant_Classification_analysis.ipynb`, and
`notebooks/analysis/OddsPath_classifications.ipynb` (and is intended for scripts that
work with the same data, not just notebooks), all of which read from and
write to a common set of directories under a caller-supplied `data_dir`
(each caller resolves its own `PROJECT_ROOT`/`data_dir` since that depends
on where it's actually running):

- `input_maves_dir` (`data/input/maves`): curated MAVE inputs, e.g.
  `CHEK2_Gebbia_2024.xlsx` and `Supplementary_Data_3.xlsx`.
- `mave_data_dir` (`data/output/maves`): the integrated variant effect
  dataset produced upstream of these three notebooks.
- `supplementary_data_dir` (`data/output/supplementary_data`):
  `Supplementary_Data_4/5/6.xlsx` -- written by
  `Variant_Classification_analysis.ipynb` (5, and reads 4) and
  `OddsPath_classifications.ipynb` (6).
- `mave_calibration_dir` (`data/output/mave_calibration`):
  `OddsPath_calibrations.csv.gz`, written by `OddsPath_calculations.ipynb`.
- `mave_calibration_oddspath_dir` (`data/output/mave_calibration/oddspath`):
  the individual per-classification calibration CSVs written by
  `OddsPath_classifications.ipynb`.
- `predictor_calibration_gene_specific_dir`
  (`data/output/predictor_calibration/gene_specific`): the individual
  gene-specific REVEL/AlphaMissense/MutPred2 calibration CSVs written by
  `Variant_Classification_analysis.ipynb`.
- `reclassification_dir` (`data/output/reclassification`): the Sankey
  dataset (`integrated_variant_effect_dataset_analysis.csv.gz`), written by
  `Variant_Classification_analysis.ipynb` and read by both it and
  `OddsPath_classifications.ipynb`.

Call `DataPaths.make_output_dirs()` once per run to create every output
directory (input directories are expected to already exist).
"""

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class DataPaths:
    input_maves_dir: Path
    mave_data_dir: Path
    supplementary_data_dir: Path
    mave_calibration_dir: Path
    mave_calibration_oddspath_dir: Path
    predictor_calibration_gene_specific_dir: Path
    reclassification_dir: Path

    def make_output_dirs(self) -> None:
        """Create every output directory (not `input_maves_dir`), idempotently."""
        for directory in (
            self.mave_data_dir,
            self.supplementary_data_dir,
            self.mave_calibration_dir,
            self.mave_calibration_oddspath_dir,
            self.predictor_calibration_gene_specific_dir,
            self.reclassification_dir,
        ):
            directory.mkdir(parents=True, exist_ok=True)


def get_data_paths(data_dir: Path) -> DataPaths:
    """Build the shared directory layout rooted at `data_dir` (a project's `data/`)."""
    mave_calibration_dir = data_dir / "output" / "mave_calibration"
    return DataPaths(
        input_maves_dir=data_dir / "input" / "maves",
        mave_data_dir=data_dir / "output" / "maves",
        supplementary_data_dir=data_dir / "output" / "supplementary_data",
        mave_calibration_dir=mave_calibration_dir,
        mave_calibration_oddspath_dir=mave_calibration_dir / "oddspath",
        predictor_calibration_gene_specific_dir=data_dir
        / "output"
        / "predictor_calibration"
        / "gene_specific",
        reclassification_dir=data_dir / "output" / "reclassification",
    )
