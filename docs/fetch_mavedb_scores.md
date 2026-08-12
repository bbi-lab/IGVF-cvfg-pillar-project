# Fetch MaveDB Scores (via the API)

`src/fetch_mavedb_scores.py` is an API-based equivalent of
`src/mavedb_scores.sql`. Both produce `cvfg_variants.0.tsv`-shaped output -- one row
per variant, with score, mapped-HGVS, and target-gene columns -- but the SQL
version queries a local mirror of MaveDB's Postgres database, while this
script hits the public API (`https://api.mavedb.org`) directly. Use it when a
DB mirror isn't available.

## Output columns

Same header as `cvfg_variants.0.tsv`: `gene_symbol`, `score_set_title`,
`variant_urn`, `variant_id`, `raw_hgvs_nt`, `raw_hgvs_pro`, `score`,
`rna_score`, `rna_score_d6`, `rna_score_d20`, `score_data`, `reference_based`,
`mavedb_mapped_hgvs_g`, `mavedb_mapped_hgvs_c`, `mavedb_mapped_hgvs_p`,
`mavedb_mapping_error`, `target_sequence_type`, `target_sequence`,
`preferred_transcript`.

## Known differences from the SQL/DB version

- **No internal integer `variant_id`.** MaveDB's public API only exposes it
  via a per-variant GET (`/api/v1/variants/{urn}`), which isn't practical to
  call once per variant across the tens of thousands of variants in this
  dataset. The column is kept, for schema parity, but left blank.
- **At most one of `mavedb_mapped_hgvs_g` / `_c` / `_p` is populated per
  row.** The DB's `mapped_variants` table stores all three simultaneously
  (computed together by MaveDB's internal mapping pipeline), but the public
  API's mapped-variant model exposes a single HGVS expression per mapping, at
  whatever alignment level MaveDB's mapper used for that score set: genomic
  (`hgvs.g`) for DNA-level score sets, protein (`hgvs.p`) for protein-level
  ones.
- The live API returns transient `504`s fairly often; the script retries
  automatically (`--retries`, default 8, exponential backoff), but a full run
  across all ~90 score sets can still take a while and may need a retry.

## Source of truth: `src/mavedb_scores.sql`

The score-set URN list and the five dataset-specific manual overrides at the
end of `src/mavedb_scores.sql` (the CHEK2 preferred transcript, the BRCA1_Findlay
MANE Select version lift, the NDUFAF6/PTEN `raw_hgvs_nt` blanking, and the
JAG1/TARDBP(x2)/SCN5A/KCNH2 target-sequence corrections) are **parsed
directly out of that SQL file at run time** rather than duplicated here as
Python literals -- several of them are thousand-character DNA/protein
sequences, so re-typing them would risk a silent transcription error.
Parsing the SQL keeps this script from drifting out of sync with it if the
overrides are ever revised. Point `--sql-source` at a different file if you
need to override this.

## Usage

```bash
poetry run python -m src.fetch_mavedb_scores \
  --output data/input/maves/cvfg_variants.0.api.tsv
```

Fetch a subset instead of the full ~90-dataset list (useful for testing
against the flaky live API), by repeating `--score-set-urn`:

```bash
poetry run python -m src.fetch_mavedb_scores \
  --output /tmp/smoke.tsv \
  --score-set-urn urn:mavedb:00000060-a-1 \
  --score-set-urn urn:mavedb:00001205-a-1
```

## CLI options

| Option | Default | Description |
|---|---|---|
| `--output` | `data/input/maves/cvfg_variants.0.api.tsv` | Where to write the resulting TSV |
| `--score-set-urn` | (full list from `src/mavedb_scores.sql`) | Fetch only this score-set URN; repeatable |
| `--sql-source` | `src/mavedb_scores.sql` | Source of the score-set list and manual overrides |
| `--retries` | `8` | HTTP retries per request on 5xx responses |
| `--timeout` | `120` | Per-request timeout, in seconds |
| `--sleep-between-requests` | `0.0` | Seconds to sleep between score sets, to go easier on the API |
