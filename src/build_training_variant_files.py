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

MutPred2: <input-dir>/mp2_actual_training_data.txt.gz
  One row per gene: column 1 is an opaque accession id, column 2 is a
  comma-separated list of 1-letter protein variants (e.g. "N697S"), column
  3 is a comma-separated list of "1"/"0" flags aligned positionally with
  column 2 (1 = that variant was actually used in MutPred2 training), column
  4 is the gene's protein sequence (unused here), and column 5 is the gene's
  NCBI/Entrez Gene ID -- *not* an HGNC ID, confirmed by cross-checking known
  genes (e.g. BRCA1 -> 672, BRCA2 -> 675, BAP1 -> 8314). There is no header
  row.

  This file covers many more genes than this project cares about, so rows
  are first filtered down to the Entrez Gene IDs of the genes curated in
  Supplementary_Data_3.xlsx's "Curation" sheet ("Gene" column -- a
  comma-separated cell for CALM1/CALM2/CALM3), via the checked-in
  <input-dir>/supplementary_data_3_gene_entrez_ids.tsv lookup (see below for
  why this can't just use Supplementary_Data_3's own "HGNC ID" column).

<input-dir>/supplementary_data_3_gene_entrez_ids.tsv
  Static gene_symbol -> entrez_gene_id lookup for every gene symbol
  currently in Supplementary_Data_3.xlsx's Curation sheet, resolved via
  mygene.info (POST to https://mygene.info/v3/query, scopes=symbol,
  species=human). Checked in rather than resolved at pipeline run time so
  this step stays offline/reproducible. Supplementary_Data_3.xlsx's own
  "HGNC ID" column can't be used for this lookup directly -- as of this
  writing it has at least 3 wrong values (DDX3X listed under CARD11's HGNC
  ID 16393; KCNQ4 listed under 6296, one off from its real 6298; SCN5A
  listed under 6331, which is actually SCN5A's *Entrez* Gene ID, not its
  HGNC ID) -- so this file is keyed on gene symbol, which is correct in
  Supplementary_Data_3.xlsx, instead.

  If Supplementary_Data_3.xlsx's Curation sheet ever gains a gene symbol not
  present here, `build_mutpred2_training_file` raises rather than silently
  dropping that gene's training variants -- regenerate this file (add the
  missing symbol's mygene.info-resolved Entrez Gene ID) when that happens.

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

import pandas as pd
from Bio.SeqUtils import seq3

_AA_CODE_RE = re.compile(r"^([A-Za-z*])(\d+)([A-Za-z*])$")

DEFAULT_INPUT_DIR = Path("data/input/predictors")
DEFAULT_OUTPUT_DIR = Path("data/intermediate/variant_annotation/data")
DEFAULT_GENE_LIST_SOURCE = Path("data/input/maves/Supplementary_Data_3.xlsx")

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


def _curated_gene_symbols(gene_list_source: Path) -> set[str]:
    curation = pd.read_excel(gene_list_source, sheet_name="Curation")
    symbols: set[str] = set()
    for cell in curation["Gene"]:
        symbols.update(symbol.strip() for symbol in str(cell).split(","))
    return symbols


def _entrez_id_to_gene_symbol(entrez_map_path: Path, gene_list_source: Path) -> dict[str, str]:
    curated_symbols = _curated_gene_symbols(gene_list_source)

    symbol_to_entrez_id: dict[str, str] = {}
    with entrez_map_path.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            symbol_to_entrez_id[row["gene_symbol"]] = row["entrez_gene_id"]

    missing = curated_symbols - symbol_to_entrez_id.keys()
    if missing:
        raise ValueError(
            f"{gene_list_source} has gene symbol(s) {sorted(missing)} with no Entrez Gene ID in "
            f"{entrez_map_path} -- regenerate that lookup (see build_training_variant_files.py's "
            "module docstring) before rerunning."
        )

    return {symbol_to_entrez_id[symbol]: symbol for symbol in curated_symbols}


def build_mutpred2_training_file(
    src_path: Path, entrez_map_path: Path, gene_list_source: Path, dest_path: Path
) -> None:
    entrez_id_to_gene_symbol = _entrez_id_to_gene_symbol(entrez_map_path, gene_list_source)

    keys: set[tuple[str, str]] = set()
    with _open_text(src_path) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            gene_symbol = entrez_id_to_gene_symbol.get(row[4].strip())
            if gene_symbol is None:
                continue
            gene_symbol = CALM_GENE_ALIASES.get(gene_symbol, gene_symbol)
            variants = row[1].split(",")
            training_flags = row[2].split(",")
            for variant, is_training in zip(variants, training_flags):
                if is_training.strip() != "1":
                    continue
                keys.add((gene_symbol, aa1_to_unqualified_hgvs_p(variant)))

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
    p.add_argument(
        "--mutpred2-source", default=None, type=Path, help="Default: <input-dir>/mp2_actual_training_data.txt.gz"
    )
    p.add_argument(
        "--mutpred2-gene-entrez-map",
        default=None,
        type=Path,
        help="Default: <input-dir>/supplementary_data_3_gene_entrez_ids.tsv",
    )
    p.add_argument(
        "--gene-list-source",
        default=DEFAULT_GENE_LIST_SOURCE,
        type=Path,
        help=f"Default: {DEFAULT_GENE_LIST_SOURCE} (Supplementary_Data_3.xlsx's 'Curation' sheet)",
    )
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
        args.mutpred2_source or input_dir / "mp2_actual_training_data.txt.gz",
        args.mutpred2_gene_entrez_map or input_dir / "supplementary_data_3_gene_entrez_ids.tsv",
        args.gene_list_source,
        args.mutpred2_dest or output_dir / "mutpred2_training_variants.tsv",
    )


if __name__ == "__main__":
    main()
