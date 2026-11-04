"""Recompute the archived AAV9 apo/holo comparison using Python's standard library."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parents[1] / 'comparison'
CONDITIONS = ('holo_complete', 'apo_complete')


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pose_models(text, names):
    """Match heavy atoms by name, retaining Vina rank and the receptor frame."""
    models, atoms, score = [], None, None
    for line in text.splitlines():
        if line.startswith('MODEL'):
            if atoms is not None or int(line.split()[1]) != len(models) + 1:
                raise ValueError('Truncated model or noncontiguous pose ranks')
            atoms, score = {}, None
        elif line.startswith('REMARK VINA RESULT:'):
            if atoms is None or score is not None:
                raise ValueError('Missing model or duplicate Vina score')
            score = float(line.split()[3])
        elif line.startswith(('ATOM  ', 'HETATM')):
            if atoms is None:
                raise ValueError('Atom outside a pose model')
            if line.split()[-1] in ('H', 'HD', 'HS'):
                continue
            name = line[12:16].strip()
            if name in atoms:
                raise ValueError('Duplicate heavy-atom name')
            atoms[name] = tuple(float(line[i:i + 8]) for i in (30, 38, 46))
        elif line.startswith('ENDMDL'):
            if atoms is None or set(atoms) != set(names) or score is None:
                raise ValueError('Incomplete heavy-atom map or Vina score')
            if not all(math.isfinite(x) for xyz in atoms.values() for x in xyz) or not math.isfinite(score):
                raise ValueError('Nonfinite pose coordinate or score')
            models.append((score, [atoms[name] for name in names]))
            atoms = None
    if atoms is not None or not models:
        raise ValueError('Empty or truncated pose file')
    return models


def rmsd(a, b):
    """Coordinate displacement in a shared receptor frame, without ligand fitting."""
    if not a or len(a) != len(b) or any(len(p) != 3 for p in [*a, *b]):
        raise ValueError('RMSD needs matching, nonempty 3D coordinate arrays')
    if not all(math.isfinite(x) for p in [*a, *b] for x in p):
        raise ValueError('Nonfinite RMSD coordinates')
    return math.sqrt(sum(math.dist(p, q) ** 2 for p, q in zip(a, b)) / len(a))


def analyze(inputs):
    sources = read_json(inputs / 'SOURCES.json')
    for name, record in sources['files'].items():
        if sha256(inputs / name) != record['sha256']:
            raise ValueError(f'Changed archived input: {name}')
    reference = read_json(inputs / 'reference.json')
    names, native = reference['names'], reference['coordinates']
    if len(names) != len(set(names)) or len(names) != len(native):
        raise ValueError('Reference atom names must form a bijection')
    rmsd(native, native)
    config, runs = read_json(inputs / 'configuration.json'), read_json(inputs / 'runs.json')
    seeds = config['seeds']
    expected = {(condition, seed) for condition in CONDITIONS for seed in seeds}
    actual = [(r['receptor'], r['seed']) for r in runs]
    if len(set(seeds)) != 3 or len(actual) != 6 or set(actual) != expected:
        raise ValueError('Comparison requires all six matched condition/seed runs')
    threshold, rows = config['primary_pose_threshold_A'], []
    for run in sorted(runs, key=lambda r: (CONDITIONS.index(r['receptor']), r['seed'])):
        if run['status'] != 'ok' or run['exit_code'] != 0:
            raise ValueError('Archived run did not complete; do not silently omit it')
        path = inputs / run['pose_file']
        if sha256(path) != run['pose_sha256']:
            raise ValueError('Pose differs from the original run record')
        models = pose_models(path.read_text(encoding='utf-8'), names)
        if len(models) > config['num_modes']:
            raise ValueError('More returned poses than the recorded configuration')
        distances = [rmsd(xyz, native) for _, xyz in models]
        rows.append({
            'condition': run['receptor'], 'seed': run['seed'], 'poses': len(models),
            'top1_score_kcal_mol': models[0][0], 'top1_displacement_A': distances[0],
            'top5_min_displacement_A': min(distances[:5]),
            'best_retained_displacement_A': min(distances),
            'retained_within_threshold': min(distances) <= threshold,
        })
    summaries = []
    for condition in CONDITIONS:
        selected = [row for row in rows if row['condition'] == condition]
        summaries.append({
            'condition': condition, 'completed_runs': len(selected), 'eligible_runs': 3,
            'retained_poses': sum(row['poses'] for row in selected),
            **{key: statistics.mean(row[key] for row in selected) for key in (
                'top1_score_kcal_mol', 'top1_displacement_A', 'top5_min_displacement_A',
                'best_retained_displacement_A')},
            'runs_retaining_pose_within_threshold': sum(row['retained_within_threshold'] for row in selected),
        })
    # Round displayed derived values only; original coordinates and scores stay intact.
    def displayed(records):
        return [{k: round(v, 3) if isinstance(v, float) else v for k, v in row.items()} for row in records]
    return {'coordinate_threshold_A': threshold, 'rows': displayed(rows),
            'summaries': displayed(summaries), 'source_commit': sources['source_commit']}


def report(result):
    header = '| State | Runs | Retained poses | Top1 score (kcal/mol) | Top1 displacement (Å) | Top5 minimum (Å) | Best retained (Å) |'
    means = [header, '| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    individual = [header.replace('Runs | Retained poses', 'Seed | Poses'), means[1]]
    def cells(row):
        return ' | '.join(f'{row[key]:.3f}' for key in (
            'top1_score_kcal_mol', 'top1_displacement_A', 'top5_min_displacement_A', 'best_retained_displacement_A'))
    for row in result['summaries']:
        means.append(f"| {row['condition'].split('_')[0]} | {row['completed_runs']}/{row['eligible_runs']} | {row['retained_poses']} | {cells(row)} |")
    for row in result['rows']:
        individual.append(f"| {row['condition'].split('_')[0]} | {row['seed']} | {row['poses']} | {cells(row)} |")
    within = sum(row['retained_within_threshold'] for row in result['rows'])
    return f"""# AAV9 apo versus holo: archived local-docking comparison

Six completed AutoDock Vina 1.2.7 runs compare ligand-stripped holo (7MTZ)
with independently solved apo (7MT0), using the same three seeds.
This is a retrospective coordinate diagnostic for one biological system.

## Results

Values below are seed means. Displacement is named-heavy-atom RMSD from
one deposited modeled GAL copy, in the aligned receptor frame.

{chr(10).join(means)}

**{within}/6 runs** retained a pose within {result['coordinate_threshold_A']:g} Å of the chosen modeled reference.
Favorable Top1 Vina scores coincide with substantial coordinate displacement
in both states. These observations do not establish experimental binding
affinity, specificity, or biologically incorrect recognition.

## Individual runs

{chr(10).join(individual)}

## What was compared

Both inputs retain three whole subunits (author chains A, B and I) around
the same reference-derived 20 Å cube. The historical `complete` suffix
means this local neighborhood; the docking input is not a full capsid.
Before cropping, apo coordinates were aligned using 31,080 corresponding Cα
atoms from the full assemblies (0.607 Å RMSD). The ligand was not fitted separately.

The archived setup used standard Vina scoring, three seeds (104729, 130363,
155921), exhaustiveness 16, two CPUs and up to 20 returned poses per run.
It used a rigid receptor and dry preparation. Input ligand coordinates were
generated independently; the search box and neighborhood used the known site.

## Interpretation and limits

The apo input has a lower mean Top1 coordinate displacement in these runs;
the holo input has a more favorable mean Top1 Vina score. Neither quantity
establishes which state models experimental recognition better.

The AAV9 modeled GAL coordinates have not passed an independent atomic-pose
reference qualification here. RMSD and the 2 Å count are descriptive coordinate
agreement, not validated pose accuracy. This local, known-site experiment
does not measure blind site discovery. Seeds are computational repeats, not
independent biological samples. One fixed reference copy is used; no capsid-wide
symmetry minimization or internal sampling diagnosis is claimed. There are
no tested nonbinders or assay measurements in this comparison.

The wider hypotheses in [PROTOCOL.md](../PROTOCOL.md) remain research proposals.
This completed comparison answers the smaller descriptive apo/holo question.

## Reproduction and sources

Run `python3 pilot/pilot-script.py --check` from the repository root.
It verifies input hashes and recomputes both this table and `results.json`
from the archived PDBQT models, using Python's standard library. It does not
rerun docking or independently reproduce receptor/ligand preparation.

[Input provenance](inputs/SOURCES.json) records each original path, SHA256
and extraction. Only the six selected runs and the reference fields needed
for this calculation were extracted; poses and logs are unchanged.
Source snapshot: `{result['source_commit']}`.
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--check', action='store_true', help='verify committed results without writing')
    mode.add_argument('--write', action='store_true', help='regenerate RESULTS.md and results.json')
    args = parser.parse_args()
    try:
        result = analyze(ROOT / 'inputs')
        outputs = {'RESULTS.md': report(result), 'results.json': json.dumps(result, indent=2, sort_keys=True) + '\n'}
        for name, text in outputs.items():
            path = ROOT / name
            if args.check and path.read_text(encoding='utf-8') != text:
                raise ValueError(f'Stale result: {path}')
            if args.write:
                path.write_text(text, encoding='utf-8')
        print('Archived inputs and both results verified.' if args.check else outputs['RESULTS.md'], end='\n' if args.check else '')
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, f'Comparison failed: {error}\n')


if __name__ == '__main__':
    main()
