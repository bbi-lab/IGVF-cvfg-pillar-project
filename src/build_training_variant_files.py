"""Regenerate revel_training_variants.tsv and mutpred2_training_variants.tsv
from richer upstream training-variant sources, for use with annotate_predictors.py's
--revel-training-file and --mutpred2-training-file.

This is a preparatory step for the variant-annotation pipeline's Step 12
(REVEL and AlphaMissense annotation; see docs/variant_annotation_pipeline.md
and docs/build_training_variant_files.md) -- run it beforehand so Step 12 can
pick up its output. `scripts/variant_annotation_pipeline.sh`'s `step_12` does
this automatically.

Reads from data/input/predictors/ (this project's own committed source data,
overridable with --input-dir) and writes into
data/intermediate/variant_annotation/data/ (the gitignored staging directory
`scripts/run_variant_annotation_pipeline.sh` mounts as the variant-annotation
pipeline's VARIANT_DATA_DIR, overridable with --output-dir; created if
missing).

Sources
-------

REVEL: <input-dir>/clustering_variants_revel_training_overlap.csv.gz
  Genome-wide REVEL-training-overlap list with genomic coordinates and a
  1-letter protein_variant column (e.g. "Q416K"). All rows are SNVs.

MutPred2: <input-dir>/MutPred2_scores_1.csv.gz and <input-dir>/MutPred2_scores_2.csv.gz
  Two independently-generated exports (Gene/AA/MP2_train and
  Gene/variant/is_mp2_train respectively) that fully agree on every
  (gene, variant) pair they share; unioned here for maximum coverage.

CALM1/CALM2/CALM3
------------------

CALM1, CALM2, and CALM3 encode an identical protein, so a MAVE assay of
"calmodulin" can't attribute a variant to one paralog over another. This
pipeline's gene_symbol column always records these variants as "CALM1"
(never a joint label), so MutPred2 training entries for CALM2/CALM3 are
aliased onto CALM1 here -- otherwise they would never match the input's
gene_symbol and the training-set overlap would silently look empty for
every CALM1 variant. REVEL matching is genomic-coordinate-based (no gene
symbol involved), so no aliasing is needed there.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path

from Bio.SeqUtils import seq3

_AA_CODE_RE = re.compile(r"^([A-Za-z*])(\d+)([A-Za-z*])$")

DEFAULT_INPUT_DIR = Path("data/input/predictors")
DEFAULT_OUTPUT_DIR = Path("data/intermediate/variant_annotation/data")

CALM_GENE_ALIASES = {"CALM2": "CALM1", "CALM3": "CALM1"}


def _open_text(path: Path):
    """Open `path` for text reading, transparently decompressing .gz sources."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="", encoding="utf-8")
    return path.open(newline="", encoding="utf-8")


def aa1_to_unqualified_hgvs_p(code: str) -> str:
    """'Q416K' -> 'p.Gln416Lys'"""
    m = _AA_CODE_RE.match(code.strip())
    if not m:
        raise ValueError(f"Unparseable 1-letter AA code: {code!r}")
    ref, pos, alt = m.groups()
    return f"p.{seq3(ref)}{pos}{seq3(alt)}"


def build_revel_training_file(src_path: Path, dest_path: Path) -> None:
    # Keyed on genomic coordinates only (gene_symbol is informational, per
    # the loader's contract), so overlapping transcripts/genes at the same
    # locus collapse to one row.
    rows: dict[tuple[str, int, int, str, str], tuple[str, str]] = {}
    with _open_text(src_path) as fh:
        for row in csv.DictReader(fh):
            chrom = row["#CHROM"]
            if chrom.startswith("chr"):
                chrom = chrom[3:]
            pos = int(row["POS"])
            key = (chrom, pos, pos, row["REF"].upper(), row["ALT"].upper())
            if key not in rows:
                rows[key] = (row["gene_symbol"], aa1_to_unqualified_hgvs_p(row["protein_variant"]))

    with dest_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["gene_symbol", "unqualified_hgvs_p", "chromosome", "hg38_start", "hg38_end", "ref_allele", "alt_allele"]
        )
        for (chrom, start, stop, ref, alt), (gene, hgvs_p) in sorted(rows.items()):
            writer.writerow([gene, hgvs_p, chrom, start, stop, ref, alt])

    print(f"Wrote {len(rows)} REVEL training-variant rows to {dest_path}")


def _add_true_rows(path: Path, gene_col: str, variant_col: str, train_col: str, into: set[tuple[str, str]]) -> None:
    with _open_text(path) as fh:
        for row in csv.DictReader(fh):
            if row[train_col].strip() != "True":
                continue
            gene = CALM_GENE_ALIASES.get(row[gene_col], row[gene_col])
            into.add((gene, aa1_to_unqualified_hgvs_p(row[variant_col])))


def build_mutpred2_training_file(scores1_path: Path, scores2_path: Path, dest_path: Path) -> None:
    keys: set[tuple[str, str]] = set()
    _add_true_rows(scores1_path, "Gene", "AA", "MP2_train", keys)
    _add_true_rows(scores2_path, "Gene", "variant", "is_mp2_train", keys)

    with dest_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_symbol", "unqualified_hgvs_p"])
        for gene, hgvs_p in sorted(keys):
            writer.writerow([gene, hgvs_p])

    print(f"Wrote {len(keys)} MutPred2 training-variant rows to {dest_path}")


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--input-dir", default=DEFAULT_INPUT_DIR, type=Path, help="Directory containing the upstream source files"
    )
    p.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        type=Path,
        help="Directory to write the training-variant TSVs into (created if missing)",
    )
    p.add_argument(
        "--revel-source",
        default=None,
        type=Path,
        help="Default: <input-dir>/clustering_variants_revel_training_overlap.csv.gz",
    )
    p.add_argument("--mutpred2-source-1", default=None, type=Path, help="Default: <input-dir>/MutPred2_scores_1.csv.gz")
    p.add_argument("--mutpred2-source-2", default=None, type=Path, help="Default: <input-dir>/MutPred2_scores_2.csv.gz")
    p.add_argument("--revel-dest", default=None, type=Path, help="Default: <output-dir>/revel_training_variants.tsv")
    p.add_argument(
        "--mutpred2-dest", default=None, type=Path, help="Default: <output-dir>/mutpred2_training_variants.tsv"
    )
    args = p.parse_args(argv)

    input_dir = args.input_dir
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    build_revel_training_file(
        args.revel_source or input_dir / "clustering_variants_revel_training_overlap.csv.gz",
        args.revel_dest or output_dir / "revel_training_variants.tsv",
    )
    build_mutpred2_training_file(
        args.mutpred2_source_1 or input_dir / "MutPred2_scores_1.csv.gz",
        args.mutpred2_source_2 or input_dir / "MutPred2_scores_2.csv.gz",
        args.mutpred2_dest or output_dir / "mutpred2_training_variants.tsv",
    )


if __name__ == "__main__":
    main()
