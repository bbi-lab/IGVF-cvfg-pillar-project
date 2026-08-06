-- Build data/input/maves/cvfg_variants.0.tsv (Step 1's input to the
-- variant-annotation pipeline -- see docs/variant_annotation_pipeline.md)
-- from a local mirror of MaveDB's Postgres database.
--
-- Run against that mirror with psql (e.g. `psql -f src/mavedb_scores.sql
-- your_mavedb_db`), then export the final SELECT's result to
-- cvfg_variants.0.tsv (e.g. via psql's `\copy (...) to '...' with (format
-- csv, delimiter E'\t', header)`).
--
-- src/fetch_mavedb_scores.py is an API-based equivalent of this script, for
-- use when a local MaveDB DB mirror isn't available -- see
-- docs/fetch_mavedb_scores.md. It parses the score-set URN list and the
-- manual overrides below directly out of this file at run time, rather than
-- duplicating them (several are thousand-character DNA/protein sequences),
-- so keep any edits here in the shapes it expects: the `ss.urn in (...)`
-- list, the `ss.urn='...'` CALM1/2/3 query, and the `update igvf_cvfg_pipeline_input set
-- ...` overrides at the end.
--
-- Structure:
--   1. Main query (below): one row per variant effect measurement, for the
--      curated list of score sets in the `ss.urn in (...)` clause, joining
--      each variant to its score set's target gene/sequence and its current
--      mapped-variant record (genomic/transcript/protein HGVS, ClinGen
--      allele ID, mapping error). Written into a temp table, igvf_cvfg_pipeline_input.
--   2. CALM1/2/3 (urn:mavedb:00000001-c-1): appended separately, since that
--      score set needs `distinct on (hgvs_pro)` handling (see its comment
--      below) that doesn't apply to any other score set here.
--   3. Manual, dataset-specific overrides: a handful of `update igvf_cvfg_pipeline_input`
--      statements correcting known MaveDB data issues (wrong preferred
--      transcript, stale transcript version, short/wrong target sequences)
--      that would otherwise cause incorrect downstream variant mapping.
--      Each is commented with why it's needed.
--   4. Final query: selects and shapes the columns that become
--      cvfg_variants.0.tsv's header.
--
-- Get all except CALM1/2/3, which requires special handling.
drop table if exists igvf_cvfg_pipeline_input;
select
  tg.mapped_hgnc_name as gene_symbol,
  ss.title as score_set_title,
  v.id as variant_id,
  v.urn as variant_urn,
  v.hgvs_nt as raw_hgvs_nt,
  v.hgvs_pro as raw_hgvs_pro,
  tg.accession_id is not null as reference_based,
  (v.data -> 'score_data' ->> 'score')::float as score,
  (v.data -> 'score_data') as score_data,
  mv.hgvs_g,
  mv.hgvs_c,
  mv.hgvs_p,
  mv.clingen_allele_id,
  mv.error_message as mapping_error,
  ts.sequence_type as target_sequence_type,
  ts.sequence as target_sequence,
  '' as preferred_transcript
into temp igvf_cvfg_pipeline_input
from
  scoresets ss
    left outer join target_genes tg on ss.id = tg.scoreset_id -- One-to-[zero or one] for these score sets
    left outer join target_sequences ts on tg.target_sequence_id = ts.id,
  variants v
    left outer join mapped_variants mv on v.id = mv.variant_id and mv."current" = true
where ss.urn in (
  -- All datasets from the originexcept CALM1/2/3
  'urn:mavedb:00000005-a-5',
  'urn:mavedb:00000005-a-6',
  'urn:mavedb:00000013-a-1',
  'urn:mavedb:00000050-a-1',
  'urn:mavedb:00000054-a-1',
  'urn:mavedb:00000060-a-1',
  'urn:mavedb:00000060-a-2',
  'urn:mavedb:00000068-0-1',
  'urn:mavedb:00000068-a-1',
  'urn:mavedb:00000068-b-1',
  'urn:mavedb:00000068-c-1',
  'urn:mavedb:00000094-a-2',
  'urn:mavedb:00000094-a-13',
  'urn:mavedb:00000094-a-3',
  'urn:mavedb:00000094-a-12',
  'urn:mavedb:00000096-a-1',
  'urn:mavedb:00000096-b-1',
  'urn:mavedb:00000097-0-2',
  'urn:mavedb:00000098-a-2',
  'urn:mavedb:00000099-a-1',
  'urn:mavedb:00000108-0-1',
  'urn:mavedb:00000108-a-1',
  'urn:mavedb:00000108-b-1',
  'urn:mavedb:00000112-a-1',
  'urn:mavedb:00000657-a-1',
  'urn:mavedb:00000657-b-1',
  'urn:mavedb:00000659-a-1',
  'urn:mavedb:00000662-0-1',
  'urn:mavedb:00000663-a-1',
  'urn:mavedb:00000665-a-1',
  'urn:mavedb:00000665-b-1',
  'urn:mavedb:00000665-c-1',
  'urn:mavedb:00000665-d-1',
  'urn:mavedb:00000673-0-1',
  'urn:mavedb:00000674-a-2',
  'urn:mavedb:00000674-b-1',
  'urn:mavedb:00000674-c-1',
  'urn:mavedb:00000675-a-1',
  'urn:mavedb:00001197-a-5',
  'urn:mavedb:00001198-a-1',
  'urn:mavedb:00001200-0-1',
  'urn:mavedb:00001200-a-1',
  'urn:mavedb:00001200-b-1',
  'urn:mavedb:00001200-c-1',
  'urn:mavedb:00001200-d-1',
  'urn:mavedb:00001200-e-1',
  'urn:mavedb:00001205-a-1',
  'urn:mavedb:00001222-a-2',
  'urn:mavedb:00001222-b-2',
  'urn:mavedb:00001224-a-1',
  'urn:mavedb:00001226-a-1',
  'urn:mavedb:00001226-b-1',
  'urn:mavedb:00001227-a-2',
  'urn:mavedb:00001228-a-1',
  'urn:mavedb:00001229-a-2',
  'urn:mavedb:00001230-a-2',
  'urn:mavedb:00001231-a-2',
  'urn:mavedb:00001233-a-1',
  'urn:mavedb:00001234-0-1',
  'urn:mavedb:00001234-a-1',
  'urn:mavedb:00001234-b-1',
  'urn:mavedb:00001234-c-1',
  'urn:mavedb:00001234-d-1',
  'urn:mavedb:00001234-e-1',
  'urn:mavedb:00001234-f-1',
  'urn:mavedb:00001234-g-1',
  'urn:mavedb:00001234-h-1',
  'urn:mavedb:00001235-a-1',
  'urn:mavedb:00001236-0-1',
  'urn:mavedb:00001242-a-1',
  'urn:mavedb:00001250-a-2',
  'urn:mavedb:00001251-a-1',
  'urn:mavedb:00001250-a-2',
  'urn:mavedb:00001254-a-1',
  'urn:mavedb:00001259-a-2',
  'urn:mavedb:00001260-a-2',
  'urn:mavedb:00001262-a-3',
  'urn:mavedb:00001263-a-2',
  'urn:mavedb:00001264-a-2',
  'urn:mavedb:00001265-a-2',
  'urn:mavedb:00001266-b-1',
  'urn:mavedb:00001267-0-2',
  'urn:mavedb:00001268-0-1',
  'urn:mavedb:00001268-a-1',
  'urn:mavedb:00001268-b-1',
  'urn:mavedb:00001268-c-1',
  -- 7 new datasets
  'urn:mavedb:00001269-a-1',
  'urn:mavedb:00001269-b-1',
  'urn:mavedb:00001269-c-1',
  'urn:mavedb:00001280-a-1',
  'urn:mavedb:00001279-a-1',
  'urn:mavedb:00001278-a-1',
  'urn:mavedb:00001278-b-1',
  'urn:mavedb:00001277-a-1'
)
and ss.id = v.scoreset_id;

-- Add CALM1/2/3, which requires a DISTINCT ON clause. This score set was uploaded at DNA-level with scores repeated, but it's a protein-level assay.
insert into igvf_cvfg_pipeline_input
select
  distinct on (v.hgvs_pro)
  tg.mapped_hgnc_name as gene_symbol,
  ss.title as score_set_title,
  v.id as variant_id,
  v.urn as variant_urn,
  null as raw_hgvs_nt,
  v.hgvs_pro as raw_hgvs_pro,
  tg.accession_id is not null as reference_based,
  (v.data -> 'score_data' ->> 'score')::float as score,
  (v.data -> 'score_data') as score_data,
  mv.hgvs_g,
  mv.hgvs_c,
  mv.hgvs_p,
  mv.clingen_allele_id,
  mv.error_message as mapping_error,
  ts.sequence_type as target_sequence_type,
  ts.sequence as target_sequence,
  '' as preferred_transcript
from
  scoresets ss
    left outer join target_genes tg on ss.id = tg.scoreset_id -- One-to-[zero or one] for these score sets
    left outer join target_sequences ts on tg.target_sequence_id = ts.id,
  variants v
    left outer join mapped_variants mv on v.id = mv.variant_id and mv."current" = true
where ss.urn = 'urn:mavedb:00000001-c-1'
and ss.id = v.scoreset_id
order by v.hgvs_pro, v.id;

-- NDUFAF6_Sung_2024 and PTEN_Mighell_2018 were uploaded to MaveDB as DNA-level score sets but should be at protein level.
-- select * from igvf_cvfg_pipeline_input v where split_part(v.variant_urn, '#', 1) in ('urn:mavedb:00000663-a-1');
update igvf_cvfg_pipeline_input
set raw_hgvs_nt = ''
where split_part(variant_urn, '#', 1) in ('urn:mavedb:00000663-a-1', 'urn:mavedb:00000054-a-1');

-- Force the mapper to choose NP_009125.1 instead of NP_001005735.1 for CHEK2.
update igvf_cvfg_pipeline_input
set preferred_transcript = 'NM_007194.4'
where split_part(variant_urn, '#', 1) = 'urn:mavedb:00001205-a-1';

-- Lift BRCA1_Findlay_2018 over to the current MANE Select transcript version. The CDS is identical, and intronic positions in the dataset are also unchanged.
update igvf_cvfg_pipeline_input
set raw_hgvs_nt = replace(raw_hgvs_nt, 'NM_007294.3', 'NM_007294.4');

-- Extend the JAG1 target sequence, since the last two NTs are in a new exon and otherwise get mapped wrongly.
update igvf_cvfg_pipeline_input
set target_sequence = 'ATGCGTTCCCCACGGACGCGCGGCCGGTCCGGGCGCCCCCTAAGCCTCCTGCTCGCCCTGCTCTGTGCCCTGCGAGCCAAGGTGTGTGGGGCCTCGGGTCAGTTCGAGTTGGAGATCCTGTCCATGCAGAACGTGAACGGGGAGCTGCAGAACGGGAACTGCTGCGGCGGCGCCCGGAACCCGGGAGACCGCAAGTGCACCCGCGACGAGTGTGACACATACTTCAAAGTGTGCCTCAAGGAGTATCAGTCCCGCGTCACGGCCGGGGGGCCCTGCAGCTTCGGCTCAGGGTCCACGCCTGTCATCGGGGGCAACACCTTCAACCTCAAGGCCAGCCGCGGCAACGACCGCAACCGCATCGTGCTGCCTTTCAGTTTCGCCTGGCCGAGGTCCTATACGTTGCTTGTGGAGGCGTGGGATTCCAGTAATGACACCGTTCAACCTGACAGTATTATTGAAAAGGCTTCTCACTCGGGCATGATCAACCCCAGCCGGCAGTGGCAGACGCTGAAGCAGAACACGGGCGTTGCCCACTTTGAGTATCAGATCCGCGTGACCTGTGATGACTACTACTATGGCTTTGGCTGCAATAAGTTCTGCCGCCCCAGAGATGACTTCTTTGGACACTATGCCTGTGACCAGAATGGCAACAAAACTTGCATGGAAGGCTGGATGGGCCCCGAATGTAACAGAGCTATTTGCCGACAAGGCTGCAGTCCTAAGCATGGGTCTTGCAAACTCCCAGGTGACTGCAGGTGCCAGTACGGCTGGCAAGGCCTGTACTGTGATAAGTGCATCCCACACCCGGGATGCGTCCACGGCATCTGTAATGAGCCCTGGCAGTGCCTCTGTGAGACCAACTGGGGCGGCCAGCTCTGTGACAAAGATCTCAATTACTGTGGGACTCATCAGCCGTGTCTCAACGGGGGAACTTGTAGCAACACAGGCCCTGACAAATATCAGTGTTCCTGCCCTGAGGGGTATTCAGGACCCAACTGTGAAATTGCTGAGCACGCCTGCCTCTCTGATCCCTGTCACAACAGAGGCAGCT'
where target_sequence = 'ATGCGTTCCCCACGGACGCGCGGCCGGTCCGGGCGCCCCCTAAGCCTCCTGCTCGCCCTGCTCTGTGCCCTGCGAGCCAAGGTGTGTGGGGCCTCGGGTCAGTTCGAGTTGGAGATCCTGTCCATGCAGAACGTGAACGGGGAGCTGCAGAACGGGAACTGCTGCGGCGGCGCCCGGAACCCGGGAGACCGCAAGTGCACCCGCGACGAGTGTGACACATACTTCAAAGTGTGCCTCAAGGAGTATCAGTCCCGCGTCACGGCCGGGGGGCCCTGCAGCTTCGGCTCAGGGTCCACGCCTGTCATCGGGGGCAACACCTTCAACCTCAAGGCCAGCCGCGGCAACGACCGCAACCGCATCGTGCTGCCTTTCAGTTTCGCCTGGCCGAGGTCCTATACGTTGCTTGTGGAGGCGTGGGATTCCAGTAATGACACCGTTCAACCTGACAGTATTATTGAAAAGGCTTCTCACTCGGGCATGATCAACCCCAGCCGGCAGTGGCAGACGCTGAAGCAGAACACGGGCGTTGCCCACTTTGAGTATCAGATCCGCGTGACCTGTGATGACTACTACTATGGCTTTGGCTGCAATAAGTTCTGCCGCCCCAGAGATGACTTCTTTGGACACTATGCCTGTGACCAGAATGGCAACAAAACTTGCATGGAAGGCTGGATGGGCCCCGAATGTAACAGAGCTATTTGCCGACAAGGCTGCAGTCCTAAGCATGGGTCTTGCAAACTCCCAGGTGACTGCAGGTGCCAGTACGGCTGGCAAGGCCTGTACTGTGATAAGTGCATCCCACACCCGGGATGCGTCCACGGCATCTGTAATGAGCCCTGGCAGTGCCTCTGTGAGACCAACTGGGGCGGCCAGCTCTGTGACAAAGATCTCAATTACTGTGGGACTCATCAGCCGTGTCTCAACGGGGGAACTTGTAGCAACACAGGCCCTGACAAATATCAGTGTTCCTGCCCTGAGGGGTATTCAGGACCCAACTGTGAAATTGCT';

-- The two TARDBP score sets cover portions of the protein, and neither starts from position 1, but the target sequence gives the whole coding region.
update igvf_cvfg_pipeline_input
set target_sequence = 'GGTAATAGCAGAGGGGGTGGAGCTGGTTTGGGAAACAATCAAGGTAGTAATATGGGTGGTGGGATGAACTTTGGTGCGTTCAGCATTAATCCAGCCATGATGGCTGCCGCCCAGGCAGCACTACAGAGCAGTTGGGGTATGATGGGCATGTTAGCCAGCCAGCAGAACCAGTCAGGCCCATCGGGTAATAACCAAAACCAAGGCAACATGCAGAGGGAGCCAAACCAGGCCTTCGGTTCTGGAAATAACTCTTATAGTGGCTCTAATTCTGGTGCAGCAATTGGTTGGGGATCAGCATCCAATGCAGGGTCGGGCAGTGGTTTTAATGGAGGCTTTGGCTCAAGCATGGATTCTAAGTCTTCTGGCTGGGGAATG'
where split_part(variant_urn, '#', 1) = 'urn:mavedb:00000060-a-1';

update igvf_cvfg_pipeline_input
set target_sequence = 'AGCAGTTGGGGTATGATGGGCATGTTAGCCAGCCAGCAGAACCAGTCAGGCCCATCGGGTAATAACCAAAACCAAGGCAACATGCAGAGGGAGCCAAACCAGGCCTTCGGTTCTGGAAATAACTCTTATAGTGGCTCTAATTCTGGTGCAGCAATTGGTTGGGGATCAGCATCCAATGCAGGGTCGGGCAGTGGTTTTAATGGAGGCTTTGGCTCAAGCATGGATTCTAAGTCTTCTGGCTGGGGAATG'
where split_part(variant_urn, '#', 1) = 'urn:mavedb:00000060-a-2';

-- The mapper picks NP_001317189.1 instead of NP_000326.2 when given the short target sequence for SCN5A.
update igvf_cvfg_pipeline_input
set target_sequence = 'YGAAVLFLLMCTFALIAHWLACIWYAIGNMEQPHMDSRIGWLHNLGDQIGKPYNSSG'
where target_sequence = 'YGAAVLFLLMCTFALIAHWL';

-- The aligner sometimes picks NP_000229.1 instead of NP_000229.1 when we give it the short target sequence for KCNH2.
update igvf_cvfg_pipeline_input
set target_sequence = 'LFRVIRLARIGRILRLIRGAKGIRTLLFALMMSLPALFNIGLLLFLVMFIYSIFGMANFAY'
where target_sequence = 'LFRVIRLARIGR';

-- Export this as cvfg_variants.0.tsv.
select
  gene_symbol,
  score_set_title,
  variant_urn,
  variant_id,
  raw_hgvs_nt,
  raw_hgvs_pro,
  score,
  coalesce((score_data ->> 'rna_score')::float, (score_data ->> 'score_rna')::float) as rna_score,
  (score_data ->> 'rna_score_d6')::float as rna_score_d6,
  (score_data ->> 'rna_score_d20')::float as rna_score_d20,
  score_data,
  reference_based,
  hgvs_g as mavedb_mapped_hgvs_g,
  hgvs_c as mavedb_mapped_hgvs_c,
  hgvs_p as mavedb_mapped_hgvs_p,
  mapping_error as mavedb_mapping_error,
  target_sequence_type,
  target_sequence,
  preferred_transcript
from igvf_cvfg_pipeline_input
order by gene_symbol, variant_id;
