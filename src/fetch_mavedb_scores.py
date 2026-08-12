#!/usr/bin/env python3
"""Fetch MaveDB score/variant data for the CVFG project via the public MaveDB REST API.

This is an API-based equivalent of `src/mavedb_scores.sql`, which
produces the same output (`cvfg_variants.0.tsv`) by querying a local mirror of
MaveDB's Postgres database. Use this script instead when a DB mirror isn't
available -- it hits https://api.mavedb.org directly.

CLI usage:

    python -m src.fetch_mavedb_scores --output data/input/maves/cvfg_variants.0.api.tsv

Pass one or more `--score-set-urn` options to fetch a subset instead of the
full ~90-dataset list (useful for testing against the flaky live API).

The manual, dataset-specific overrides applied at the end of `mavedb_scores.sql`
(target sequence corrections for JAG1/TARDBP/SCN5A/KCNH2, the CHEK2 preferred
transcript, the BRCA1_Findlay MANE Select version lift, and the NDUFAF6/PTEN
`raw_hgvs_nt` blanking) are parsed directly out of that SQL file at run time,
rather than duplicated here as Python literals. Several of them are
thousand-character DNA/protein sequences, so re-typing them would risk a
silent transcription error; parsing the SQL keeps this script and the SQL
version from drifting apart if the overrides are ever revised.

Known differences from the SQL/DB version:

- No internal integer `variant_id`: MaveDB's public API only exposes it via
  a per-variant GET, which isn't practical to call once per variant across
  tens of thousands of variants. The output column is kept (for schema
  parity with `cvfg_variants.0.tsv`) but left blank.
- At most one of `mavedb_mapped_hgvs_g` / `_c` / `_p` is populated per row.
  The DB's `mapped_variants` table stores all three simultaneously, but the
  public API's mapped-variant model exposes a single HGVS expression per
  mapping, at whatever alignment level MaveDB's mapper used for that score
  set (genomic for DNA-level score sets, protein for protein-level ones).
"""

import io
import json
import re
import time
from pathlib import Path

import click
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

API_BASE = "https://api.mavedb.org/api/v1"
DEFAULT_SQL_SOURCE = Path(__file__).resolve().parent / "mavedb_scores.sql"
DEFAULT_OUTPUT = Path("data/input/maves/cvfg_variants.0.api.tsv")

FIXED_SCORE_COLUMNS = {"variant_urn", "raw_hgvs_nt", "hgvs_splice", "raw_hgvs_pro"}

OUTPUT_COLUMNS = [
    "gene_symbol",
    "score_set_title",
    "variant_urn",
    "variant_id",
    "raw_hgvs_nt",
    "raw_hgvs_pro",
    "score",
    "rna_score",
    "rna_score_d6",
    "rna_score_d20",
    "score_data",
    "reference_based",
    "mavedb_mapped_hgvs_g",
    "mavedb_mapped_hgvs_c",
    "mavedb_mapped_hgvs_p",
    "mavedb_mapping_error",
    "target_sequence_type",
    "target_sequence",
    "preferred_transcript",
]

# These tolerate optional whitespace around `=`/parens/commas (`\s*`) so that
# reformatting mavedb_scores.sql's whitespace/alignment doesn't silently break
# parsing; `\s+` marks separators that must remain (e.g. between keywords).
_URN_IN_CLAUSE_RE = re.compile(r"ss\.urn\s+in\s*\((.*?)\)\s+and\s+ss\.id\s*=\s*v\.scoreset_id", re.DOTALL)
_URN_LITERAL_RE = re.compile(r"'(urn:mavedb:[^']+)'")
_CALM_URN_RE = re.compile(r"ss\.urn\s*=\s*'(urn:mavedb:[^']+)'")
_BLANK_NT_RE = re.compile(
    r"set\s+raw_hgvs_nt\s*=\s*''\s+where\s+split_part\(variant_urn,\s*'#',\s*1\)\s+in\s*\((.*?)\)\s*;", re.DOTALL
)
_PREFERRED_TRANSCRIPT_RE = re.compile(
    r"set\s+preferred_transcript\s*=\s*'([^']+)'\s+where\s+split_part\(variant_urn,\s*'#',\s*1\)\s*=\s*'([^']+)'"
)
_HGVS_NT_REPLACE_RE = re.compile(r"set\s+raw_hgvs_nt\s*=\s*replace\(raw_hgvs_nt,\s*'([^']+)',\s*'([^']+)'\)")
_TARGET_SEQUENCE_BLOCK_RE = re.compile(
    r"set\s+target_sequence\s*=\s*'([A-Za-z]+)'\s+"
    r"where\s+(?:target_sequence\s*=\s*'([A-Za-z]+)'|split_part\(variant_urn,\s*'#',\s*1\)\s*=\s*'(urn:mavedb:[^']+)')"
)


def parse_score_set_urns(sql_text):
    """Extract the main `ss.urn in (...)` score-set URN list from `mavedb_scores.sql`.

    This excludes the CALM1/2/3 score set, which is handled by a separate
    query in the SQL (see `parse_calm_urn`).
    """
    match = _URN_IN_CLAUSE_RE.search(sql_text)
    if not match:
        raise ValueError("Could not find the `ss.urn in (...)` score-set list in the SQL source.")
    return _URN_LITERAL_RE.findall(match.group(1))


def parse_calm_urn(sql_text):
    """Extract the CALM1/2/3 score-set URN from `mavedb_scores.sql`'s second query."""
    match = _CALM_URN_RE.search(sql_text)
    if not match:
        raise ValueError("Could not find the CALM1/2/3 `ss.urn='...'` query in the SQL source.")
    return match.group(1)


def parse_manual_overrides(sql_text):
    """Parse the dataset-specific manual overrides applied at the end of `mavedb_scores.sql`.

    Returns a dict with keys:
      - blank_raw_hgvs_nt_urns: list of score-set URNs whose raw_hgvs_nt should be blanked
      - preferred_transcript: (transcript, urn) tuple
      - hgvs_nt_replacement: (old_substring, new_substring) tuple
      - target_sequence_replacements: {old_sequence: new_sequence} keyed by exact sequence match
      - target_sequence_by_urn: {score_set_urn: new_sequence}
    """
    blank_nt_match = _BLANK_NT_RE.search(sql_text)
    if not blank_nt_match:
        raise ValueError("Could not find the raw_hgvs_nt blanking override in the SQL source.")
    blank_raw_hgvs_nt_urns = _URN_LITERAL_RE.findall(blank_nt_match.group(1))

    preferred_transcript_match = _PREFERRED_TRANSCRIPT_RE.search(sql_text)
    if not preferred_transcript_match:
        raise ValueError("Could not find the preferred_transcript override in the SQL source.")
    preferred_transcript = preferred_transcript_match.groups()

    hgvs_nt_replace_match = _HGVS_NT_REPLACE_RE.search(sql_text)
    if not hgvs_nt_replace_match:
        raise ValueError("Could not find the raw_hgvs_nt version-lift override in the SQL source.")
    hgvs_nt_replacement = hgvs_nt_replace_match.groups()

    target_sequence_replacements = {}
    target_sequence_by_urn = {}
    for match in _TARGET_SEQUENCE_BLOCK_RE.finditer(sql_text):
        new_sequence, old_sequence, urn = match.groups()
        if old_sequence is not None:
            target_sequence_replacements[old_sequence] = new_sequence
        else:
            target_sequence_by_urn[urn] = new_sequence
    if not target_sequence_replacements and not target_sequence_by_urn:
        raise ValueError("Could not find any target_sequence overrides in the SQL source.")

    return {
        "blank_raw_hgvs_nt_urns": blank_raw_hgvs_nt_urns,
        "preferred_transcript": preferred_transcript,
        "hgvs_nt_replacement": hgvs_nt_replacement,
        "target_sequence_replacements": target_sequence_replacements,
        "target_sequence_by_urn": target_sequence_by_urn,
    }


def apply_manual_overrides(df, overrides):
    """Apply the parsed manual overrides to a combined variant dataframe, mirroring the
    `update igvf_cvfg_pipeline_input ...` statements at the end of `mavedb_scores.sql`."""
    df = df.copy()
    base_urn = df["variant_urn"].str.split("#").str[0]

    blank_nt_urns = set(overrides["blank_raw_hgvs_nt_urns"])
    df.loc[base_urn.isin(blank_nt_urns), "raw_hgvs_nt"] = ""

    transcript, pt_urn = overrides["preferred_transcript"]
    df.loc[base_urn == pt_urn, "preferred_transcript"] = transcript

    old_nt, new_nt = overrides["hgvs_nt_replacement"]
    df["raw_hgvs_nt"] = df["raw_hgvs_nt"].str.replace(old_nt, new_nt, regex=False)

    for old_sequence, new_sequence in overrides["target_sequence_replacements"].items():
        df.loc[df["target_sequence"] == old_sequence, "target_sequence"] = new_sequence

    for urn, new_sequence in overrides["target_sequence_by_urn"].items():
        df.loc[base_urn == urn, "target_sequence"] = new_sequence

    return df


def build_session(retries=8, backoff_factor=2.0):
    """Build a `requests.Session` that retries the transient 5xx errors this API returns often."""
    session = requests.Session()
    retry = Retry(
        total=retries,
        backoff_factor=backoff_factor,
        status_forcelist=(500, 502, 503, 504),
        allowed_methods=("GET",),
        respect_retry_after_header=True,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.headers.update({"User-Agent": "igvf-cvfg-pillar-project/fetch_mavedb_scores"})
    return session


def _get(session, url, *, params=None, timeout=120):
    response = session.get(url, params=params, timeout=timeout)
    response.raise_for_status()
    return response


def fetch_score_set_metadata(session, urn, *, timeout=120):
    """Fetch a score set's title and (first) target gene metadata."""
    data = _get(session, f"{API_BASE}/score-sets/{urn}", timeout=timeout).json()
    target_genes = data.get("targetGenes") or []
    target_gene = target_genes[0] if target_genes else {}
    target_sequence = target_gene.get("targetSequence") or {}
    return {
        "gene_symbol": target_gene.get("mappedHgncName") or target_gene.get("name"),
        "score_set_title": data.get("title"),
        "reference_based": target_gene.get("targetAccession") is not None,
        "target_sequence_type": target_sequence.get("sequenceType"),
        "target_sequence": target_sequence.get("sequence"),
    }


def _jsonable(value):
    return value.item() if hasattr(value, "item") else value


def fetch_scores(session, urn, *, timeout=120):
    """Fetch a score set's per-variant scores as a dataframe, one row per variant.

    Mirrors `(v.data->'score_data'->>'score')::float as score` plus the raw
    `score_data` jsonb blob from the SQL version: every score column beyond
    the fixed accession/HGVS ones is packed into a `score_data` JSON string,
    and `rna_score`/`rna_score_d6`/`rna_score_d20` are pulled out of it the
    same way the SQL does.
    """
    response = _get(session, f"{API_BASE}/score-sets/{urn}/scores", timeout=timeout)
    return scores_from_csv_text(response.text, urn=urn)


def scores_from_csv_text(csv_text, *, urn="<score set>"):
    """Parse a `/score-sets/{urn}/scores` CSV response into the scores dataframe.

    Split out from `fetch_scores` so the CSV-to-dataframe logic can be
    exercised without a live HTTP call.
    """
    df = pd.read_csv(
        io.StringIO(csv_text),
        dtype={"accession": str, "hgvs_nt": str, "hgvs_splice": str, "hgvs_pro": str},
        na_values=["NA"],
    )
    df = df.rename(columns={"accession": "variant_urn", "hgvs_nt": "raw_hgvs_nt", "hgvs_pro": "raw_hgvs_pro"})
    if "score" not in df.columns:
        raise click.ClickException(f"{urn}: /scores response has no 'score' column.")

    score_columns = [c for c in df.columns if c not in FIXED_SCORE_COLUMNS]
    df["score_data"] = [
        json.dumps({k: _jsonable(v) for k, v in row.dropna().items()}) for _, row in df[score_columns].iterrows()
    ]

    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df["rna_score"] = _first_present_numeric(df, "rna_score", "score_rna")
    df["rna_score_d6"] = _first_present_numeric(df, "rna_score_d6")
    df["rna_score_d20"] = _first_present_numeric(df, "rna_score_d20")
    return df


def _first_present_numeric(df, *column_names):
    for name in column_names:
        if name in df.columns:
            return pd.to_numeric(df[name], errors="coerce")
    return pd.Series(float("nan"), index=df.index)


def fetch_mapped_variants(session, urn, *, timeout=180):
    """Fetch a score set's current mapped variants, extracting the mapped HGVS
    expression(s), ClinGen allele ID, and mapping error per variant.

    Only `current=True` records are kept, matching the SQL's
    `mv."current"=TRUE` join condition.
    """
    records = _get(session, f"{API_BASE}/score-sets/{urn}/mapped-variants", timeout=timeout).json()
    return mapped_variants_from_records(records)


def mapped_variants_from_records(records):
    """Parse `/score-sets/{urn}/mapped-variants` JSON records into the mapped-variants dataframe.

    Split out from `fetch_mapped_variants` so the JSON-to-dataframe logic can
    be exercised without a live HTTP call.
    """
    rows = []
    for record in records:
        if not record.get("current"):
            continue
        expressions = ((record.get("postMapped") or {}).get("expressions")) or []
        hgvs_by_syntax = {expr["syntax"]: expr["value"] for expr in expressions if expr.get("syntax")}
        rows.append(
            {
                "variant_urn": record.get("variantUrn"),
                "mavedb_mapped_hgvs_g": hgvs_by_syntax.get("hgvs.g"),
                "mavedb_mapped_hgvs_c": hgvs_by_syntax.get("hgvs.c"),
                "mavedb_mapped_hgvs_p": hgvs_by_syntax.get("hgvs.p"),
                "clingen_allele_id": record.get("clingenAlleleId"),
                "mavedb_mapping_error": record.get("errorMessage"),
            }
        )
    columns = [
        "variant_urn",
        "mavedb_mapped_hgvs_g",
        "mavedb_mapped_hgvs_c",
        "mavedb_mapped_hgvs_p",
        "clingen_allele_id",
        "mavedb_mapping_error",
    ]
    return pd.DataFrame(rows, columns=columns)


def fetch_score_set(session, urn, *, is_calm=False, timeout=120):
    """Fetch and assemble one score set's rows in the `cvfg_variants.0.tsv` shape."""
    metadata = fetch_score_set_metadata(session, urn, timeout=timeout)
    scores = fetch_scores(session, urn, timeout=timeout)
    mapped = fetch_mapped_variants(session, urn, timeout=timeout)

    df = scores.merge(mapped, on="variant_urn", how="left")
    for key, value in metadata.items():
        df[key] = value
    df["variant_id"] = ""
    df["preferred_transcript"] = ""

    if is_calm:
        # This score set was uploaded at DNA level with scores repeated, but it's a
        # protein-level assay: force raw_hgvs_nt blank and dedupe on hgvs_pro, keeping
        # the lowest-numbered variant per distinct protein change (mirrors the SQL's
        # `distinct on (v.hgvs_pro) ... order by v.hgvs_pro, v.id`).
        df["raw_hgvs_nt"] = None
        variant_number = df["variant_urn"].str.rsplit("#", n=1).str[-1].astype(int)
        df = (
            df.assign(_variant_number=variant_number)
            .sort_values(["raw_hgvs_pro", "_variant_number"])
            .drop_duplicates(subset="raw_hgvs_pro", keep="first")
            .drop(columns="_variant_number")
        )

    return df


def fetch_all_score_sets(session, urns, calm_urn=None, *, timeout=120, show_progress=True, sleep_between=0.0):
    """Fetch every score set in `urns` (plus `calm_urn`, if given) and concatenate the results."""
    tasks = [(urn, False) for urn in urns]
    if calm_urn is not None:
        tasks.append((calm_urn, True))

    frames = []

    def _run(pairs):
        for urn, is_calm in pairs:
            frames.append(fetch_score_set(session, urn, is_calm=is_calm, timeout=timeout))
            if sleep_between:
                time.sleep(sleep_between)

    if show_progress:
        with click.progressbar(
            tasks, label="Fetching MaveDB score sets", item_show_func=lambda t: t[0] if t else ""
        ) as bar:
            _run(bar)
    else:
        _run(tasks)

    return pd.concat(frames, ignore_index=True)


def finalize_output(df):
    """Apply final column selection/formatting to match `cvfg_variants.0.tsv`'s shape."""
    df = df.copy()
    df["reference_based"] = df["reference_based"].map({True: "true", False: "false"})
    df = df.sort_values(["gene_symbol", "variant_urn"])
    return df[OUTPUT_COLUMNS]


@click.command()
@click.option(
    "--output",
    type=click.Path(dir_okay=False, path_type=Path),
    default=DEFAULT_OUTPUT,
    show_default=True,
    help="Where to write the resulting TSV.",
)
@click.option(
    "--score-set-urn",
    "score_set_urns",
    multiple=True,
    help="Fetch only this score-set URN. Repeatable. Defaults to the full list from mavedb_scores.sql.",
)
@click.option(
    "--sql-source",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=DEFAULT_SQL_SOURCE,
    show_default=True,
    help="mavedb_scores.sql, the source of truth for the score-set list and manual overrides.",
)
@click.option("--retries", default=8, show_default=True, help="HTTP retries per request on 5xx responses.")
@click.option("--timeout", default=120, show_default=True, help="Per-request timeout, in seconds.")
@click.option(
    "--sleep-between-requests",
    default=0.0,
    show_default=True,
    help="Seconds to sleep between score sets, to go easier on the API.",
)
def main(output, score_set_urns, sql_source, retries, timeout, sleep_between_requests):
    """Fetch MaveDB score/variant data for the CVFG project, via the MaveDB API."""
    sql_text = sql_source.read_text()
    overrides = parse_manual_overrides(sql_text)
    calm_urn = parse_calm_urn(sql_text)

    if score_set_urns:
        urns = [urn for urn in score_set_urns if urn != calm_urn]
        include_calm = calm_urn in score_set_urns
    else:
        urns = parse_score_set_urns(sql_text)
        include_calm = True

    session = build_session(retries=retries)
    combined = fetch_all_score_sets(
        session,
        urns,
        calm_urn if include_calm else None,
        timeout=timeout,
        sleep_between=sleep_between_requests,
    )
    combined = apply_manual_overrides(combined, overrides)
    result = finalize_output(combined)

    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, sep="\t", index=False, na_rep="")
    click.echo(f"Wrote {len(result)} rows to {output}")


if __name__ == "__main__":
    main()
