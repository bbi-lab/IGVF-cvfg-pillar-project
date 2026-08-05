# Recalculate ClinGen Classification

`src/recalculate_clingen_classification.py` adds two columns to a CVFG
variants TSV, recomputing each row's ClinGen Evidence Repository (erepo)
classification after discarding functional-assay evidence codes:

| Column | Description |
|---|---|
| `Updated_Classification_ClinGen_repo` | ACMG/AMP classification recomputed from the remaining evidence codes |
| `Updated_Evidence Codes_ClinGen_repo` | The remaining evidence codes themselves (comma-separated) -- what the classification was computed from |

"Functional-assay evidence codes" means PS3 and BS3 (functional studies) and
PP3 and BP4 (computational/predictor evidence), including their
`_Supporting`/`_Moderate`/`_Strong`/`_Very_Strong` strength variants -- so the
recomputed classification reflects clinical/population evidence only.

## Source logic

Both columns are ported from
`notebooks/analysis/Integrated_variant_effect_dataset_pipeline.ipynb`'s
`filter_and_recalculate` (which strips the codes and re-joins what's left)
and `classify_acmg` (which combines evidence-code strengths into an ACMG/AMP
call) functions, applied to the erepo-annotated
`clingen_evidence_repository.Applied Evidence Codes (Met)` column produced by
`annotate_erepo.py` (step 8 of the variant-annotation pipeline; see
`docs/variant_annotation_pipeline.md`).

That notebook worked from a flat erepo dump -- one row per classification,
evidence codes as a single comma-separated string. `annotate_erepo.py`
instead produces a **pipe-delimited** column, aligned one segment per
`mapped_hgvs_c` DNA candidate (the same convention `flag_variants.py`'s
`Flag` column uses). This script applies the ported logic independently to
each candidate's segment and re-joins the two results with `|` in the same
order. When a segment holds evidence from more than one matching erepo
record (`annotate_erepo.py` merges those with `" | "` when, e.g., two expert
panels classified the same variant), the codes from all of that candidate's
records are pooled before filtering and reclassifying.

A candidate with no erepo match (empty segment) is treated as having no
evidence and gets `VUS` / an empty evidence string, same as
`filter_and_recalculate` on the notebook's `NaN` case.

### A ported quirk: BA1 never counts as benign-standalone

The notebook's `classify_acmg` checks evidence strength against the literal
string `"Benign_standalone"` (lowercase `s`), but its own evidence-code table
maps `BA1` to `"Benign_Standalone"` (capital `S`). The comparison therefore
never matches, so `BA1` alone never triggers the benign-standalone
classification rule. This script preserves that behavior rather than fixing
it, so its output matches the classifications already published from that
notebook -- see `classify_acmg`'s docstring in the script.

## Usage

Locally (with the Poetry environment):

```bash
poetry run python -m src.recalculate_clingen_classification \
  data/cvfg_variants.16.tsv data/cvfg_variants.16-reclassified.tsv
```

Via Docker (same image as `flag_variants`, see `compose.yaml`):

```bash
src/scripts/run_recalculate_clingen_classification.sh \
  data/cvfg_variants.16.tsv data/cvfg_variants.16-reclassified.tsv
```

Like `run_flag_variants.sh`, this wrapper maps its input/output paths against
the `/work` staging mount (`${VARIANT_DATA_DIR:-./data}`), since it reads
from the same staged pipeline data those files live in rather than a
committed repo file.

## CLI options

| Option | Default | Description |
|---|---|---|
| `--evidence-codes-col` | `clingen_evidence_repository.Applied Evidence Codes (Met)` | Column with pipe-delimited erepo "Applied Evidence Codes (Met)" values |

Raises `click.ClickException` if `--evidence-codes-col` isn't present in the
input file.
