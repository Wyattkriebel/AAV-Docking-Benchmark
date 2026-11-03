"""Read-only replay of source-pinned sequence/domain overlap; no docking inputs.

The default verifier needs only Python's standard library. Full optimal-alignment
replay additionally requires the separately pinned BioPython 1.86 environment.
Sequence screen hits are warnings, not automatic evolutionary homology claims.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import re
from pathlib import Path

ARCHIVE = Path('data/archive/family_audit')
LEDGER = Path('benchmarks/family_overlap_audit.json')
STANDARD_AA = frozenset('ACDEFGHIKLMNPQRSTVWY')


def _load(path):
    path = Path(path)
    raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix == '.gz' else raw)


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def parse_training_rows(text: str) -> list[dict]:
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith(('#', 'PARAM.', 'HEADER:')):
            continue
        fields = line.split()
        if len(fields) != 2:
            raise ValueError(f'Invalid dataset row: {line}')
        match = re.search(r'_([0-9][A-Za-z0-9]{3})([A-Za-z0-9]+)\.pdb$', fields[1])
        if not match:
            raise ValueError(f'Missing PDB/chain identity: {fields[1]}')
        rows.append(dict(pdb_id=match[1].upper(), author_chain_hint=match[2],
                         protein_path=fields[1], prediction_path=fields[0]))
    return rows


def alignment_metrics(left: str, right: str, coordinates: list[list[int]]) -> dict:
    """Count identities and coverage from monotone alignment breakpoint arrays."""
    if len(coordinates) != 2 or len(coordinates[0]) != len(coordinates[1]):
        raise ValueError('Invalid alignment coordinates')
    if not coordinates[0]:
        return dict(aligned_nongap_pairs=0, identities=0, identity=0., shorter_coverage=0.,
                    alignment_columns_including_gaps=0, identity_including_gap_columns=0.)
    if not left or not right:
        raise ValueError('Empty sequence')
    for axis, sequence in zip(coordinates, [left, right]):
        if any(not isinstance(n, int) or n < 0 or n > len(sequence) for n in axis):
            raise ValueError('Out-of-range alignment coordinate')
    pairs = identities = columns = 0
    for i in range(len(coordinates[0]) - 1):
        l0, l1 = coordinates[0][i:i + 2]
        r0, r1 = coordinates[1][i:i + 2]
        dl, dr = l1 - l0, r1 - r0
        if dl < 0 or dr < 0 or dl == dr == 0 or (dl and dr and dl != dr):
            raise ValueError('Nonmonotone or non-diagonal alignment block')
        columns += max(dl, dr)
        if dl and dr:
            pairs += dl
            identities += sum(left[a] == right[b] and left[a] in STANDARD_AA
                              for a, b in zip(range(l0, l1), range(r0, r1)))
    return dict(aligned_nongap_pairs=pairs, identities=identities,
                identity=identities / pairs if pairs else 0.,
                shorter_coverage=pairs / min(len(left), len(right)),
                alignment_columns_including_gaps=columns,
                identity_including_gap_columns=identities / columns if columns else 0.)


def screen_hit(metrics: dict, policy: dict) -> bool:
    rule = policy['sequence_screen']
    return (metrics['identity'] >= rule['identity_min_inclusive'] and
            metrics['shorter_coverage'] >= rule['shorter_sequence_coverage_min_inclusive'])


def _aligner(policy):
    import Bio
    from Bio import Align
    from Bio.Align import substitution_matrices
    if Bio.__version__ != '1.86':
        raise ValueError('Full replay requires pinned BioPython 1.86')
    rule = policy['sequence_screen']
    aligner = Align.PairwiseAligner(mode='local')
    aligner.substitution_matrix = substitution_matrices.load(rule['substitution_matrix'])
    aligner.open_gap_score = rule['gap_open']
    aligner.extend_gap_score = rule['gap_extend']
    aligner.epsilon = rule['epsilon']
    return aligner


def _normalise(sequence, alphabet):
    return ''.join(a if a in alphabet else 'X' for a in sequence)


def _classification_matches(query, entities, scope, key):
    wanted = set(query[key])
    return [dict(record_id=e['record_id'], pdb_id=e['pdb_id'], description=e['description'],
                 matching_ids=sorted(wanted & set(e[key])),
                 training_entity_matches_listed_author_hint=e['training_entity_matches_listed_author_hint'])
            for e in entities if e['scope'] == scope and e['pdb_id'] != query['pdb_id']
            and wanted & set(e[key])]


def _warning_summary(pdb_id, scope, alignments, records):
    found = []
    for i, item in enumerate(alignments):
        if not item['homology_screen_hit']:
            continue
        left, right = records[item['left_id']], records[item['right_id']]
        if left['pdb_id'] == pdb_id:
            query, target = left, right
        elif right['pdb_id'] == pdb_id:
            query, target = right, left
        else:
            continue
        if target['scope'] != scope:
            continue
        metrics = alignment_metrics(left['sequence'], right['sequence'], item['coordinates'])
        found.append(dict(alignment_index=i, query_id=query['record_id'], target_id=target['record_id'],
                          target_parent_id=target['parent_record_id'], query_length=len(query['sequence']),
                          target_length=len(target['sequence']), identity=item['identity'],
                          shorter_coverage=item['shorter_coverage'], aligned_nongap_pairs=item['aligned_nongap_pairs'],
                          identity_including_gap_columns=metrics['identity_including_gap_columns'],
                          score=item['score'],
                          training_entity_matches_listed_author_hint=target['training_entity_matches_listed_author_hint']))
    found.sort(key=lambda x: (-x['identity'], -x['aligned_nongap_pairs'], x['query_id'], x['target_id']))
    return dict(warning_alignment_count=len(found),
                warning_target_parent_ids=sorted({x['target_parent_id'] for x in found}),
                top_identity_warnings=found[:10], all_warning_alignment_indices=sorted(x['alignment_index'] for x in found),
                interpretation='Prespecified percent/coverage screen warnings; short, gap-rich, repeat-driven and cross-fold matches are not automatically credible evolutionary homologs.')


def summarize(root: Path) -> dict:
    archive = root / ARCHIVE
    policy = _load(archive / 'qualification_policy.json')
    entities = _load(archive / 'protein_entities.json.gz')
    records = {r['record_id']: r for r in _load(archive / 'sequence_records.json.gz')}
    alignments = _load(archive / 'all_alignments.json.gz')
    retrieval = _load(archive / 'retrieval_manifest.json')
    run = _load(archive / 'alignment_run.json')
    missing = {scope: [e['record_id'] for e in entities if e['scope'] == scope and not e['cath_superfamilies']]
               for scope in ['candidate', 'development', 'p2rank_training']}
    cases = []
    for query in entities:
        if query['scope'] != 'candidate' or query['entity_id'] != '1':
            continue
        comparisons = {}
        for scope in ['development', 'p2rank_training', 'candidate']:
            warning = _warning_summary(query['pdb_id'], scope, alignments, records)
            cath = _classification_matches(query, entities, scope, 'cath_superfamilies')
            topology = _classification_matches(query, entities, scope, 'cath_topologies')
            query_scop = {x['id'] for x in query['scop2_superfamily_lineage_nodes']}
            scop = [dict(record_id=e['record_id'], matching_nodes=[x for x in e['scop2_superfamily_lineage_nodes'] if x['id'] in query_scop])
                    for e in entities if e['scope'] == scope and e['pdb_id'] != query['pdb_id']
                    and query_scop & {x['id'] for x in e['scop2_superfamily_lineage_nodes']}]
            if cath:
                status = 'seen_CATH_homologous_superfamily'
            elif warning['warning_alignment_count']:
                status = 'unknown_sequence_screen_warning'
            elif not query['cath_superfamilies']:
                status = 'unknown_missing_candidate_CATH'
            elif missing[scope]:
                status = 'no_detected_overlap_unknown_comparison_CATH_coverage'
            else:
                status = 'novel_under_declared_screen'
            comparisons[scope] = dict(status=status, CATH_homologous_superfamily_matches=cath,
                                      CATH_topology_matches=topology, SCOP2_superfamily_matches=scop,
                                      sequence_warning_summary=warning)
        cases.append(dict(pdb_id=query['pdb_id'], candidate_entity=query['record_id'],
                          canonical_sequence_sha256=query['canonical_sequence_sha256'], length=query['length'],
                          CATH_homologous_superfamilies=query['cath_superfamilies'],
                          SCOP2_superfamilies=query['scop2_superfamily_lineage_nodes'],
                          domain_boundaries_source='data/archive/family_audit/protein_entities.json.gz#' + query['record_id'],
                          strict_novel_family_certified=False, learned_method_training_heldout_certified=False,
                          comparisons=comparisons))
    return dict(schema_version=1, source_only=True, prediction_outcomes_consulted=False,
                status='executed_source_sequence_and_domain_overlap_audit; no universal novelty certification',
                policy_archive=str(ARCHIVE / 'qualification_policy.json'), policy_sha256=_sha(archive / 'qualification_policy.json'),
                retrieval_manifest=str(ARCHIVE / 'retrieval_manifest.json'),
                retrieval_manifest_sha256=_sha(archive / 'retrieval_manifest.json'),
                candidate_inventory_sha256=retrieval['candidate_inventory_sha256'],
                alignment_software=run, training_dataset=retrieval['dataset'],
                exact_candidate_training_PDB_overlap=retrieval['exact_candidate_training_pdb_overlap'],
                exact_development_training_PDB_overlap=retrieval['exact_development_training_pdb_overlap'],
                source_classification_coverage=dict(
                    protein_entity_counts={scope: sum(e['scope'] == scope for e in entities) for scope in missing},
                    missing_CATH_entity_ids=missing,
                    selected_training_entities_without_CATH=[e['record_id'] for e in entities if e['scope'] == 'p2rank_training'
                                                             and e['training_entity_matches_listed_author_hint'] and not e['cath_superfamilies']],
                    interpretation='A missing CATH assignment or boundary is unknown. Negative sequence screening does not supply the missing structural classification.'),
                summary_counts=dict(candidate_count=len(cases),
                                    development_CATH_overlap=sum(bool(c['comparisons']['development']['CATH_homologous_superfamily_matches']) for c in cases),
                                    development_sequence_warning=sum(bool(c['comparisons']['development']['sequence_warning_summary']['warning_alignment_count']) for c in cases),
                                    P2Rank_training_CATH_overlap=sum(bool(c['comparisons']['p2rank_training']['CATH_homologous_superfamily_matches']) for c in cases),
                                    strict_novel_families_certified=0, learned_method_training_heldouts_certified=0),
                limitations=[
                    'CATH annotations are the versions integrated in the pinned RCSB snapshots, predominantly v4_3_0; absence does not imply a novel fold.',
                    'All deposited protein entities and source domains were screened; training-PDB co-entities can overflag exposure. Each hit records whether its entity matches the dataset author-chain hint.',
                    'Percent identity uses nongap residue-residue columns as frozen before calculation. Gap-inclusive identity is reported descriptively; it is not an alternative outcome-dependent acceptance rule.',
                    'Full chains and source-defined CATH/SCOP2 domains include non-recognition domains. Recognition-domain assignment and ambiguous sequence warnings need independent source adjudication before a family split is certified.',
                    'CHEN11-fpocket is 251 rows and 241 unique PDB entries, not 11 proteins. This audit does not audit score-calibration/feature-table provenance, GNINA or all historical learned-model training.',
                    'ProtGlycanDock exact train overlap for 3ZWE and valid overlap for 1I3H remain as recorded in the candidate inventory; this P2Rank audit cannot clear them.'],
                candidates=cases)


def verify(root: Path, recompute_alignments: bool = False) -> dict:
    archive = root / ARCHIVE
    pins = _load(archive / 'artifact_manifest.json')
    for item in pins['files']:
        if _sha(root / item['archive']) != item['sha256']:
            raise ValueError('Family-audit checksum mismatch: ' + item['archive'])
    retrieval = _load(archive / 'retrieval_manifest.json')
    if _sha(archive / 'qualification_policy.json') != retrieval['policy_sha256']:
        raise ValueError('Policy changed after source retrieval')
    if _sha(root / 'benchmarks/heldout_candidate_inventory.json') != retrieval['candidate_inventory_sha256']:
        raise ValueError('Candidate inventory changed')
    if _sha(archive / 'chen11-fpocket.ds') != retrieval['dataset']['sha256']:
        raise ValueError('Official training dataset changed')
    rows = parse_training_rows((archive / 'chen11-fpocket.ds').read_text())
    if rows != _load(archive / 'dataset_rows.json') or len(rows) != 251 or len({r['pdb_id'] for r in rows}) != 241:
        raise ValueError('Training row/PDB denominator mismatch')
    source_entities = {}
    source_ids = set()
    for batch in retrieval['batches']:
        if _sha(root / batch['request_archive']) != batch['request_sha256']:
            raise ValueError('RCSB request changed')
        raw_gz = (root / batch['response_archive']).read_bytes()
        raw = gzip.decompress(raw_gz)
        if hashlib.sha256(raw_gz).hexdigest() != batch['response_compressed_sha256'] or hashlib.sha256(raw).hexdigest() != batch['response_sha256']:
            raise ValueError('RCSB response changed')
        data = json.loads(raw)
        if data.get('errors'):
            raise ValueError('Unresolved GraphQL errors')
        entries = data['data']['entries']
        if {e['rcsb_id'] for e in entries} != set(batch['pdb_ids']):
            raise ValueError('Missing RCSB entry')
        source_ids.update(e['rcsb_id'] for e in entries)
        for entry in entries:
            for entity in entry['polymer_entities'] or []:
                if entity['entity_poly']['rcsb_entity_polymer_type'] == 'Protein':
                    source_entities[entity['rcsb_id']] = entity
    entities = _load(archive / 'protein_entities.json.gz')
    if set(source_entities) != {e['record_id'] for e in entities}:
        raise ValueError('Protein entity inventory differs from sources')
    for entity in entities:
        source = source_entities[entity['record_id']]
        sequence = ''.join(source['entity_poly']['pdbx_seq_one_letter_code_can'].split())
        cath = sorted({a['annotation_id'] for i in source['polymer_entity_instances'] or []
                       for a in i['rcsb_polymer_instance_annotation'] or [] if a['type'] == 'CATH'
                       and re.fullmatch(r'\d+\.\d+\.\d+\.\d+', a['annotation_id'])})
        if sequence != entity['canonical_sequence'] or cath != entity['cath_superfamilies']:
            raise ValueError('Sequence/CATH source derivation mismatch')
    records = {r['record_id']: r for r in _load(archive / 'sequence_records.json.gz')}
    entity_map = {e['record_id']: e for e in entities}
    for record in records.values():
        parent = entity_map[record['parent_record_id']]
        if record['kind'] == 'full_chain':
            expected = parent['canonical_sequence']
        else:
            positions = sorted({n for start, end in record['segments'] for n in range(start, end + 1)
                                if 1 <= n <= parent['length']})
            expected = ''.join(parent['canonical_sequence'][n - 1] for n in positions)
            domains = parent['domain_annotations']
            if not any(d['type'] == record['source_domain_type'] and d['classification_id'] == record['classification_id']
                       and d['segments'] == record['segments'] for d in domains):
                raise ValueError('Domain record lacks pinned boundary evidence')
        if expected != record['sequence'] or hashlib.sha256(expected.encode()).hexdigest() != record['sequence_sha256']:
            raise ValueError('Domain/full-chain sequence mismatch')
    policy = _load(archive / 'qualification_policy.json')
    queries = [r for r in records.values() if r['scope'] == 'candidate' and r['entity_id'] == '1']
    targets = [r for r in records.values() if r['scope'] != 'candidate'] + queries
    expected_pairs = {tuple(sorted([q['record_id'], t['record_id']])) for q in queries for t in targets if q['pdb_id'] != t['pdb_id']}
    alignments = _load(archive / 'all_alignments.json.gz')
    if len(alignments) != len(expected_pairs) or {(a['left_id'], a['right_id']) for a in alignments} != expected_pairs:
        raise ValueError('Pairwise alignment denominator mismatch')
    aligner = _aligner(policy) if recompute_alignments else None
    for item in alignments:
        left, right = records[item['left_id']]['sequence'], records[item['right_id']]['sequence']
        metrics = alignment_metrics(left, right, item['coordinates'])
        for key in ['aligned_nongap_pairs', 'identities', 'identity', 'shorter_coverage']:
            if abs(metrics[key] - item[key]) > 1e-12:
                raise ValueError('Alignment identity/coverage mismatch')
        if screen_hit(metrics, policy) != item['homology_screen_hit'] or not math.isfinite(item['score']):
            raise ValueError('Alignment screen classification mismatch')
        if aligner is not None:
            alpha = set(aligner.substitution_matrix.alphabet)
            optimal = aligner.align(_normalise(left, alpha), _normalise(right, alpha))
            coords = optimal[0].coordinates.tolist() if optimal.score > 0 else [[], []]
            if float(optimal.score) != item['score'] or coords != item['coordinates']:
                raise ValueError('Optimal affine-gap alignment replay mismatch')
    expected_summary = summarize(root)
    if expected_summary != _load(root / LEDGER):
        raise ValueError('Family overlap summary is stale')
    return dict(status='verified', source_only=True, prediction_outcomes_consulted=False,
                source_archive_count=len(pins['files']), RCSB_PDB_entries=len(source_ids),
                protein_entities=len(entities), sequence_records=len(records), pairwise_alignments=len(alignments),
                all_alignment_metrics_recomputed=True, optimal_alignments_recomputed=recompute_alignments,
                summary_counts=expected_summary['summary_counts'])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--recompute-alignments', action='store_true')
    args = parser.parse_args(argv)
    try:
        report = verify(args.root, args.recompute_alignments)
    except (OSError, KeyError, ValueError, ImportError) as exc:
        print(json.dumps(dict(status='failed', error=str(exc))))
        return 1
    print(json.dumps(report, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
