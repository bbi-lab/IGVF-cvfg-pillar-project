import json

import pandas as pd
import pytest

from src.fetch_mavedb_scores import (
    apply_manual_overrides,
    finalize_output,
    mapped_variants_from_records,
    parse_calm_urn,
    parse_manual_overrides,
    parse_score_set_urns,
    scores_from_csv_text,
)

# A small stand-in for mavedb_scores.sql's structure, not its real content -- keeps
# these tests independent of future edits to the real overrides.
FIXTURE_SQL = """
drop table if exists igvf_cvfg_pipeline_input;
select
  tg.mapped_hgnc_name gene_symbol
into temp igvf_cvfg_pipeline_input
from
  scoresets ss
where ss.urn in
(
    'urn:mavedb:00000005-a-5', 'urn:mavedb:00000005-a-6',
    'urn:mavedb:00000013-a-1'
)
and ss.id=v.scoreset_id
;

insert into igvf_cvfg_pipeline_input
select tg.mapped_hgnc_name gene_symbol
from scoresets ss
where ss.urn='urn:mavedb:00000001-c-1'
and ss.id=v.scoreset_id
;

update igvf_cvfg_pipeline_input set raw_hgvs_nt='' where split_part(variant_urn, '#', 1) in ('urn:mavedb:00000663-a-1', 'urn:mavedb:00000054-a-1');

update igvf_cvfg_pipeline_input
set preferred_transcript='NM_007194.4'
where split_part(variant_urn, '#', 1)='urn:mavedb:00001205-a-1';

update igvf_cvfg_pipeline_input
set raw_hgvs_nt=replace(raw_hgvs_nt, 'NM_007294.3', 'NM_007294.4');

update igvf_cvfg_pipeline_input
set target_sequence='ATGNEWSEQ'
where target_sequence='ATGOLDSEQ';

update igvf_cvfg_pipeline_input
set target_sequence='ATGTARDBPNEW'
where split_part(variant_urn, '#', 1)='urn:mavedb:00000060-a-1';
"""


def test_parse_score_set_urns_excludes_calm():
    urns = parse_score_set_urns(FIXTURE_SQL)
    assert urns == ["urn:mavedb:00000005-a-5", "urn:mavedb:00000005-a-6", "urn:mavedb:00000013-a-1"]


def test_parse_calm_urn():
    assert parse_calm_urn(FIXTURE_SQL) == "urn:mavedb:00000001-c-1"


def test_parse_manual_overrides():
    overrides = parse_manual_overrides(FIXTURE_SQL)
    assert overrides["blank_raw_hgvs_nt_urns"] == ["urn:mavedb:00000663-a-1", "urn:mavedb:00000054-a-1"]
    assert overrides["preferred_transcript"] == ("NM_007194.4", "urn:mavedb:00001205-a-1")
    assert overrides["hgvs_nt_replacement"] == ("NM_007294.3", "NM_007294.4")
    assert overrides["target_sequence_replacements"] == {"ATGOLDSEQ": "ATGNEWSEQ"}
    assert overrides["target_sequence_by_urn"] == {"urn:mavedb:00000060-a-1": "ATGTARDBPNEW"}


def test_parse_manual_overrides_missing_section_raises():
    with pytest.raises(ValueError, match="preferred_transcript"):
        parse_manual_overrides(FIXTURE_SQL.replace("set preferred_transcript", "set something_else"))


def test_apply_manual_overrides():
    overrides = parse_manual_overrides(FIXTURE_SQL)
    df = pd.DataFrame(
        {
            "variant_urn": [
                "urn:mavedb:00000663-a-1#1",
                "urn:mavedb:00001205-a-1#1",
                "urn:mavedb:00000060-a-1#1",
                "urn:mavedb:00000097-0-2#1",
            ],
            "raw_hgvs_nt": ["NM_001.1:c.1A>T", "", "", "NM_007294.3:c.5565A>T"],
            "target_sequence": ["ATGXYZ", "ATGXYZ", "ATGOLD_UNRELATED", "ATGOLDSEQ"],
            "preferred_transcript": ["", "", "", ""],
        }
    )

    result = apply_manual_overrides(df, overrides)

    # NDUFAF6/PTEN-style blanking, scoped to its own URN only.
    assert result.loc[0, "raw_hgvs_nt"] == ""
    # CHEK2 preferred-transcript override, scoped to its own URN only.
    assert result.loc[1, "preferred_transcript"] == "NM_007194.4"
    assert result.loc[0, "preferred_transcript"] == ""
    # TARDBP target_sequence override, scoped by URN regardless of its prior sequence value.
    assert result.loc[2, "target_sequence"] == "ATGTARDBPNEW"
    # BRCA1_Findlay version lift, global substring replace.
    assert result.loc[3, "raw_hgvs_nt"] == "NM_007294.4:c.5565A>T"
    # JAG1-style target_sequence override, matched by exact old-sequence value.
    assert result.loc[3, "target_sequence"] == "ATGNEWSEQ"


def test_apply_manual_overrides_target_sequence_by_urn_wins_over_unrelated_rows():
    overrides = parse_manual_overrides(FIXTURE_SQL)
    df = pd.DataFrame(
        {
            "variant_urn": ["urn:mavedb:00000060-a-1#1", "urn:mavedb:00000060-a-1#2"],
            "raw_hgvs_nt": ["", ""],
            "target_sequence": ["whatever-it-was-before", "whatever-it-was-before"],
            "preferred_transcript": ["", ""],
        }
    )
    result = apply_manual_overrides(df, overrides)
    assert (result["target_sequence"] == "ATGTARDBPNEW").all()


def test_scores_from_csv_text_builds_score_data_and_flattened_columns():
    csv_text = (
        "accession,hgvs_nt,hgvs_splice,hgvs_pro,score,rna_score_d6\n"
        "urn:mavedb:00000001-a-1#1,NA,NA,p.Gly1Phe,0.5,1.5\n"
        "urn:mavedb:00000001-a-1#2,NA,NA,p.Gly1Trp,-0.25,NA\n"
    )

    df = scores_from_csv_text(csv_text, urn="urn:mavedb:00000001-a-1")

    assert list(df["variant_urn"]) == ["urn:mavedb:00000001-a-1#1", "urn:mavedb:00000001-a-1#2"]
    assert df["score"].tolist() == [0.5, -0.25]
    assert df["rna_score_d6"].tolist()[0] == 1.5
    assert pd.isna(df["rna_score_d6"].tolist()[1])
    assert pd.isna(df["rna_score"]).all()  # no rna_score/score_rna column in this fixture

    first_score_data = json.loads(df.loc[0, "score_data"])
    assert first_score_data == {"score": 0.5, "rna_score_d6": 1.5}


def test_scores_from_csv_text_coalesces_score_rna_alias():
    csv_text = "accession,hgvs_nt,hgvs_splice,hgvs_pro,score,score_rna\nurn:mavedb:x#1,NA,NA,p.Gly1Phe,0.1,2.2\n"
    df = scores_from_csv_text(csv_text)
    assert df["rna_score"].tolist() == [2.2]


def test_scores_from_csv_text_requires_score_column():
    from click import ClickException

    with pytest.raises(ClickException):
        scores_from_csv_text("accession,hgvs_nt,hgvs_splice,hgvs_pro\nurn:mavedb:x#1,NA,NA,p.Gly1Phe\n")


def test_mapped_variants_from_records_filters_current_and_extracts_hgvs():
    records = [
        {
            "variantUrn": "urn:mavedb:x#1",
            "current": False,
            "postMapped": {"expressions": [{"syntax": "hgvs.g", "value": "NC_1.1:g.1A>T"}]},
            "clingenAlleleId": "CA1",
            "errorMessage": None,
        },
        {
            "variantUrn": "urn:mavedb:x#1",
            "current": True,
            "postMapped": {"expressions": [{"syntax": "hgvs.p", "value": "NP_1.1:p.Gly1Phe"}]},
            "clingenAlleleId": "CA2",
            "errorMessage": None,
        },
        {
            "variantUrn": "urn:mavedb:x#2",
            "current": True,
            "postMapped": None,
            "clingenAlleleId": None,
            "errorMessage": "mapping failed",
        },
    ]

    df = mapped_variants_from_records(records)

    assert len(df) == 2
    row1 = df[df["variant_urn"] == "urn:mavedb:x#1"].iloc[0]
    assert row1["mavedb_mapped_hgvs_p"] == "NP_1.1:p.Gly1Phe"
    assert row1["mavedb_mapped_hgvs_g"] is None
    assert row1["clingen_allele_id"] == "CA2"

    row2 = df[df["variant_urn"] == "urn:mavedb:x#2"].iloc[0]
    assert row2["mavedb_mapping_error"] == "mapping failed"


def test_mapped_variants_from_records_empty():
    df = mapped_variants_from_records([])
    assert list(df.columns) == [
        "variant_urn",
        "mavedb_mapped_hgvs_g",
        "mavedb_mapped_hgvs_c",
        "mavedb_mapped_hgvs_p",
        "clingen_allele_id",
        "mavedb_mapping_error",
    ]
    assert len(df) == 0


def test_finalize_output_orders_columns_and_stringifies_reference_based():
    df = pd.DataFrame(
        {
            "gene_symbol": ["BBB", "AAA"],
            "score_set_title": ["t2", "t1"],
            "variant_urn": ["urn:mavedb:b#1", "urn:mavedb:a#1"],
            "variant_id": ["", ""],
            "raw_hgvs_nt": ["", ""],
            "raw_hgvs_pro": ["p.1", "p.2"],
            "score": [1.0, 2.0],
            "rna_score": [None, None],
            "rna_score_d6": [None, None],
            "rna_score_d20": [None, None],
            "score_data": ["{}", "{}"],
            "reference_based": [True, False],
            "mavedb_mapped_hgvs_g": [None, None],
            "mavedb_mapped_hgvs_c": [None, None],
            "mavedb_mapped_hgvs_p": [None, None],
            "mavedb_mapping_error": [None, None],
            "target_sequence_type": ["dna", "dna"],
            "target_sequence": ["ATG", "ATG"],
            "preferred_transcript": ["", ""],
        }
    )

    result = finalize_output(df)

    assert list(result["gene_symbol"]) == ["AAA", "BBB"]  # sorted by gene_symbol, variant_urn
    assert list(result["reference_based"]) == ["false", "true"]
    assert result.columns[0] == "gene_symbol"
    assert result.columns[-1] == "preferred_transcript"
