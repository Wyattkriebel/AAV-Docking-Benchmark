"""Analytic controls for the source-only family overlap screen."""
from pathlib import Path

import pytest

from aav_benchmark.family_audit import alignment_metrics, parse_training_rows, screen_hit, verify


POLICY = {'sequence_screen': {'identity_min_inclusive': .30,
                              'shorter_sequence_coverage_min_inclusive': .70}}


def test_dataset_parser_preserves_chain_and_counts_unique_pdbs():
    text = '''# benchmark
PARAM.PREDICTION_METHOD=fpocket
HEADER: prediction protein
out1.pdb chen11/a.001.001.001_1s69a.pdb
out2.pdb chen11/a.001.001.002_1s69b.pdb
out3.pdb chen11/a.001.001.003_2hyrA.pdb
'''
    rows = parse_training_rows(text)
    assert [(r['pdb_id'], r['author_chain_hint']) for r in rows] == [('1S69', 'a'), ('1S69', 'b'), ('2HYR', 'A')]
    assert len(rows) == 3 and len({r['pdb_id'] for r in rows}) == 2


def test_dataset_parser_rejects_unidentified_source():
    with pytest.raises(ValueError, match='Missing PDB'):
        parse_training_rows('out.pdb arbitrary.pdb')


def test_identity_and_coverage_do_not_count_insertions_as_pairs():
    result = alignment_metrics('ACDEFGHIKL', 'ACDEYYFGHIKL', [[0, 4, 4, 10], [0, 4, 6, 12]])
    assert result['identities'] == result['aligned_nongap_pairs'] == 10
    assert result['identity'] == result['shorter_coverage'] == 1
    assert result['alignment_columns_including_gaps'] == 12
    assert result['identity_including_gap_columns'] == pytest.approx(10 / 12)


def test_gapped_shorter_span_is_not_coverage():
    result = alignment_metrics('ACDEFGH', 'ACDEYYYYYY', [[0, 4, 7], [0, 4, 4]])
    assert result['aligned_nongap_pairs'] == 4
    assert result['shorter_coverage'] == pytest.approx(4 / 7)
    assert not screen_hit(result, POLICY)


def test_ambiguous_residues_never_supply_identities():
    result = alignment_metrics('AXBZ*', 'AXBZ*', [[0, 5], [0, 5]])
    assert result['identities'] == 1
    assert result['identity'] == .2
    assert not screen_hit(result, POLICY)


def test_screen_thresholds_are_inclusive_and_joint():
    assert screen_hit({'identity': .3, 'shorter_coverage': .7}, POLICY)
    assert not screen_hit({'identity': .299999, 'shorter_coverage': 1}, POLICY)
    assert not screen_hit({'identity': 1, 'shorter_coverage': .699999}, POLICY)


@pytest.mark.parametrize('coordinates', [[[0, 11], [0, 10]], [[0, 4, 3], [0, 4, 5]], [[0, 4], [0, 3]]])
def test_invalid_alignment_coordinates_are_rejected(coordinates):
    with pytest.raises(ValueError):
        alignment_metrics('ACDEFGHIKL', 'ACDEFGHIKL', coordinates)


def test_archived_audit_recomputes_complete_denominator_and_decisions():
    root = Path(__file__).resolve().parents[1]
    report = verify(root)
    assert report['RCSB_PDB_entries'] == 262
    assert report['protein_entities'] == 297
    assert report['sequence_records'] == 1365
    assert report['pairwise_alignments'] == 56366
    assert report['summary_counts']['development_CATH_overlap'] == 2
    assert report['summary_counts']['P2Rank_training_CATH_overlap'] == 7
    assert report['summary_counts']['strict_novel_families_certified'] == 0
    assert report['prediction_outcomes_consulted'] is False
