"""Pinned source-only paired HADDOCK adapter and foreground execution driver.

Preparation reads frozen box metadata, independent ligands and audited protein
topologies only. No accuracy/reference coordinates are loaded. Native execution
requires an external root freeze matching the prepared plan and this source.
"""
from __future__ import annotations

import argparse
import datetime
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

import numpy as np

from .haddock_cohort_preparation import ROOT, pin, _checked
from .haddock_method import pdb_atoms, read_cns_topology

SOURCE_MANIFEST = Path('results/structural_pose_v3/vina_carb_prepared_inputs.json')
SOURCE_SHA = '0dc9e3975d5451d3e21b2f474b7218ecf4d82e110e197fab95630bb13c22579b'
PREFLIGHT_RECEIPT = Path('results/haddock_prospective_preparation/preparation_receipt.json')
PREFLIGHT_SHA = '9acaf59029d153b66764bacb8ba485ad26b7c27cc097b57237e1e2d8f62756cf'
PROPOSAL = Path('benchmarks/haddock_paired_cohort_v2_proposal.json')
ASSEMBLY_POLICY = Path('benchmarks/haddock_assembly_protocol_amendment.json')
DEFAULT_OUTPUT = Path('results/haddock_paired_cohort')
SEEDS = (104729, 130363, 155921)


class PairedCohortError(ValueError):
    pass


def _write(path, data):
    Path(path).write_text(json.dumps(data, indent=2) + '\n')


def source_box(condition, chemistry, method):
    """Use original top1 candidate even when the Vina-specific prep failed."""
    if method == 'oracle':
        jobs = condition.get('jobs', [])
        if not jobs:
            return None
        center = jobs[0]['center_A']
    else:
        candidates = condition.get('candidate_selection', {}).get('candidates', [])
        if not candidates:
            return None
        candidate = next((c for c in candidates if c.get('rank') == 1), None)
        if candidate is None:
            raise PairedCohortError('Missing original rank-one candidate')
        center = candidate.get('center_A', candidate.get('center'))
    size = chemistry['box_size']
    if np.isscalar(size):
        size = [size] * 3
    if len(center) != 3 or len(size) != 3 or not all(math.isfinite(float(x)) for x in list(center) + list(size)) or min(size) <= 0:
        raise PairedCohortError('Invalid frozen box')
    return dict(center=list(center), size=list(size))


def box_air(protein_atoms, glycan_atoms, topology, box):
    """A single atom-union AIR; no residue expansion or native contact input."""
    center, half = np.array(box['center']), np.array(box['size']) / 2
    selected = [a for a in protein_atoms if a['element'] not in ['H', 'D'] and
                np.all(np.abs(np.array(a['xyz_A']) - center) <= half)]
    if not selected:
        raise PairedCohortError('empty_box_atom_selection')
    typed = {(a['segment'], a['resnum'], a['resname'], a['atom_name']) for a in topology['atoms'].values() if a['mass'] > 2.5}
    def identity(atom):
        key = tuple(atom[k] for k in ['chain', 'resnum', 'resname', 'atom_name'])
        if key not in typed:
            raise PairedCohortError('Protein AIR identity absent from verified CNS topology')
        return key
    protein_keys = sorted(identity(a) for a in selected)
    ligand_keys = sorted(tuple(a[k] for k in ['chain', 'resnum', 'resname', 'atom_name']) for a in glycan_atoms if a['element'] not in ['H', 'D'])
    if not ligand_keys or any(k[0] != 'L' for k in ligand_keys):
        raise PairedCohortError('Incomplete or incorrectly segmented glycan selection')
    protein_selections = [(k[0], k[1], k[3]) for k in protein_keys]
    ligand_selections = [(k[0], k[1], k[3]) for k in ligand_keys]
    if len(protein_keys) != len(set(protein_selections)) or len(ligand_keys) != len(set(ligand_selections)):
        raise PairedCohortError('Ambiguous CNS atom identity')
    def selection(keys):
        lines = []
        for chain, number, _, name in keys:
            if not re.fullmatch(r'[A-Za-z0-9]', chain) or not re.fullmatch(r'-?\d+', number) or not re.fullmatch(r'[A-Za-z0-9]+', name):
                raise PairedCohortError('Unsupported CNS selection identifier')
            lines.append(f'  (segid {chain} and resid {number} and name {name})')
        return '\n or\n'.join(lines)
    text = f"assign (\n{selection(protein_keys)}\n) (\n{selection(ligand_keys)}\n) 2.0 2.0 0.0\n"
    return text, dict(box=box, protein_heavy_atom_keys=[list(k) for k in protein_keys],
                      full_glycan_heavy_atom_keys=[list(k) for k in ligand_keys], assignment_count=1,
                      native_contacts_used=False, residue_expansion=False, random_removal=False,
                      interpretation='Coarse ambiguous patch-contact pseudo-distance; not a hard Cartesian box or exact Vina search-space equivalence.')


def configuration(seed, charged_cter):
    """Native modules only; all initial/refined models retained without clustering."""
    terminal = '\n[topoaa.mol1]\ncharged_cter = true\n' if charged_cter else ''
    return ('run_dir = "run"\nmode = "local"\nncores = 1\ndebug = true\nclean = false\n'
            'molecules = ["../data/receptor.pdb", "../data/glycan.pdb"]\n'
            '[topoaa]\niniseed = 917\nautotoppar = false\n' + terminal +
            f'\n[rigidbody]\niniseed = {seed}\nsampling = 20\ntolerance = 20\nw_vdw = 1\n'
            'ambig_fname = "../data/box_air.tbl"\nrandremoval = false\n'
            '\n[flexref]\nsampling_factor = 1\ntolerance = 20\n'
            'ambig_fname = "../data/box_air.tbl"\nrandremoval = false\n'
            'mdsteps_rigid = 500\nmdsteps_cool1 = 500\nmdsteps_cool2 = 1000\nmdsteps_cool3 = 1000\n')


def assembly_protocol_support(protein_atoms):
    """Default native modules treat each segment as an independent body."""
    chains = sorted({a['chain'] for a in protein_atoms if a['element'] not in ['H', 'D']})
    return dict(source_receptor_chains=chains, supported=len(chains) == 1,
                criterion='exactly_one_source_receptor_chain_for_this_bounded_native_default_adapter',
                interpretation='Adapter/workflow assembly-preservation limitation, not a general HADDOCK multi-chain capability claim.')


def prepare_plan(root=ROOT, output=None):
    root = Path(root); output = Path(output or root / DEFAULT_OUTPUT)
    if (output / 'prepared_plan.json').exists():
        raise PairedCohortError('Refuse to overwrite prepared cohort plan')
    _checked(root / SOURCE_MANIFEST, SOURCE_SHA)
    _checked(root / PREFLIGHT_RECEIPT, PREFLIGHT_SHA)
    source = json.loads((root / SOURCE_MANIFEST).read_text())
    preflight = json.loads((root / PREFLIGHT_RECEIPT).read_text())
    proposal = json.loads((root / PROPOSAL).read_text())
    by_profile = {(r['profile'], r['case_id']): r for r in preflight['case_profiles']}
    conditions, jobs = [], []
    for row in source['conditions']:
        profile, case, method = row['receptor_profile'], row['case_id'], row['method']
        condition_id = row['condition_id']; destination = output / profile / case / method
        destination.mkdir(parents=True, exist_ok=True)
        source_pin = row['source_condition']; _checked(root / source_pin['path'], source_pin['sha256'])
        condition = json.loads((root / source_pin['path']).read_text())
        profile_dir = root / 'data/processed/heldout_pose_v3' / profile / case
        chemistry_path = profile_dir / 'chemistry.json'
        chemistry = json.loads(chemistry_path.read_text())
        assembly = assembly_protocol_support(pdb_atoms(profile_dir / 'receptor.pdb'))
        box = source_box(condition, chemistry, method)
        support = by_profile[profile, case]
        record = dict(condition_id=condition_id, profile=profile, case_id=case, method=method,
                      eligible_workflow_condition=True, source_condition=source_pin,
                      chemistry=pin(chemistry_path, root), box=box, source_preflight=support,
                      assembly_protocol=assembly,
                      information='native-derived frozen box center only' if method == 'oracle' else 'frozen protein-only predictor top1 box',
                      status='ready')
        if box is None:
            record['status'] = 'predictor_unavailable' if method != 'oracle' else 'oracle_box_unavailable'
        elif not support['glycan_supported'] or not support['exact_receptor_supported']:
            record['status'] = 'source_parameterization_gate_failed'
        elif not assembly['supported']:
            record['status'] = 'assembly_protocol_unsupported'
        else:
            audit_path = root / support['audit']['path']; _checked(audit_path, support['audit']['sha256'])
            preflight_dir = audit_path.parent
            data = destination / 'data'; data.mkdir(exist_ok=True)
            for name in ['receptor.pdb', 'glycan.pdb', 'glycan_translated.sdf']:
                shutil.copyfile(preflight_dir / 'data' / name, data / name)
            shutil.copyfile(preflight_dir / 'preparation.json', data / 'preflight_preparation.json')
            protein = pdb_atoms(data / 'receptor.pdb'); glycan = pdb_atoms(data / 'glycan.pdb')
            topology_path = preflight_dir / 'run-topology/0_topoaa/receptor_haddock.psf'
            topology = read_cns_topology(topology_path)
            try:
                air, air_metadata = box_air(protein, glycan, topology, box)
            except PairedCohortError as error:
                if str(error) != 'empty_box_atom_selection':
                    raise
                record['status'] = 'empty_box_atom_selection'
            else:
                (data / 'box_air.tbl').write_text(air)
                air_metadata.update(source_condition=source_pin, receptor=pin(data / 'receptor.pdb', root),
                                    verified_CNS_topology=pin(topology_path, root), AIR=pin(data / 'box_air.tbl', root))
                _write(data / 'box_air.json', air_metadata)
                record['AIR'] = pin(data / 'box_air.tbl', root)
                record['AIR_metadata'] = pin(data / 'box_air.json', root)
                record['input_files'] = [pin(data / name, root) for name in ['receptor.pdb', 'glycan.pdb', 'glycan_translated.sdf', 'preflight_preparation.json', 'box_air.tbl', 'box_air.json']]
        conditions.append(record)
        for seed in SEEDS:
            job = dict(job_id=condition_id + '/s' + str(seed), condition_id=condition_id, seed=seed,
                       profile=profile, case_id=case, method=method, eligible_workflow_pool=True,
                       status=record['status'], initial_model_slots=20, maximum_refined_model_slots=20,
                       model_seed_range=[seed + 1, seed + 20])
            if record['status'] == 'ready':
                directory = destination / ('s' + str(seed)); directory.mkdir(exist_ok=True)
                (directory / 'run.cfg').write_text(configuration(seed, support['charged_cter']))
                job.update(directory=str(directory.relative_to(root)), config=pin(directory / 'run.cfg', root), input_files=record['input_files'])
            jobs.append(job)
    if len(conditions) != 36 or len(jobs) != 108 or len({j['job_id'] for j in jobs}) != 108:
        raise PairedCohortError('Fixed36-condition/108-cell denominator mismatch')
    lock = json.loads((root / 'benchmarks/haddock.lock.json').read_text())
    runtime_root = root / 'tools/haddock3-2026.10.1/env/lib/python3.11/site-packages/haddock'
    runtime_sources = [pin(p, root) for p in sorted(runtime_root.rglob('*'))
                       if p.is_file() and p.suffix in {'.py', '.cns', '.top', '.param', '.yaml', '.cfg'}]
    plan = dict(schema_version=1, status='prepared_source_only_waiting_for_root_freeze', source_only_preparation=True,
                native_cohort_docking_started=False, no_accuracy_evaluation=True,
                proposal=pin(root / PROPOSAL, root), driver=pin(Path(__file__), root),
                assembly_protocol_policy=pin(root / ASSEMBLY_POLICY, root),
                source_manifest=pin(root / SOURCE_MANIFEST, root), preflight_receipt=pin(root / PREFLIGHT_RECEIPT, root),
                release_lock=pin(root / 'benchmarks/haddock.lock.json', root), CNS_binary=lock['cns_binary'],
                release_distribution=lock['distribution'], ncores=1, seeds=list(SEEDS),
                runtime_sources=runtime_sources, native_environment_freeze=lock['environment_freeze'],
                condition_denominator=36, seed_cell_denominator=108, maximum_initial_models=2160,
                maximum_refined_models=2160, native_modules=['topoaa', 'rigidbody', 'flexref'],
                ranking='Final per-seed individual HADDOCK score ascending, tie immutable originalmodelnumber; no cluster selection or native reference.',
                conditions=conditions, jobs=jobs)
    _write(output / 'prepared_plan.json', plan)
    return plan


def _pdb_text(path):
    return gzip.open(path, 'rt').read() if path.suffix == '.gz' else path.read_text()


def model_inventory(run_dir, job_id):
    """Export method scores and immutable model paths, never pose accuracy."""
    stages = []
    for directory in sorted(Path(run_dir).glob('[0-9]*_*')):
        if not directory.is_dir():
            continue
        io = directory / 'io.json'
        outputs = json.loads(io.read_text()).get('output', []) if io.is_file() else []
        by_name = {o.get('file_name'): o for o in outputs if isinstance(o, dict) and o.get('file_name')}
        rows = []
        for pdb in sorted(set(directory.glob('*.pdb')) | set(directory.glob('*.pdb.gz'))):
            name = pdb.name.removesuffix('.gz'); output = by_name.get(name, {})
            number = int(re.search(r'_(\d+)\.pdb$', name).group(1)) if re.search(r'_(\d+)\.pdb$', name) else None
            rows.append(dict(model_id=job_id + '/' + directory.name + '/' + name, original_model_number=number,
                             file_name=pdb.name, sha256=hashlib.sha256(pdb.read_bytes()).hexdigest(),
                             atom_count=sum(l.startswith(('ATOM  ', 'HETATM')) for l in _pdb_text(pdb).splitlines()),
                             score=output.get('score'), energies=output.get('unw_energies'), actual_CNS_seed=output.get('seed'),
                             cluster_id=output.get('clt_id'), cluster_rank=output.get('clt_rank'), within_cluster_rank=output.get('clt_model_rank')))
        ranked = sorted([r for r in rows if r['score'] is not None and math.isfinite(r['score'])], key=lambda r: (r['score'], r['original_model_number'] or 0, r['model_id']))
        for rank, row in enumerate(ranked, 1):
            row['individual_score_rank'] = rank
        stages.append(dict(stage=directory.name, actual_physical_model_files=len(rows), models=rows,
                           io_sha256=hashlib.sha256(io.read_bytes()).hexdigest() if io.is_file() else None,
                           individual_score_order=[r['model_id'] for r in ranked]))
    return stages


def verify_execution_freeze(plan_path, freeze_path, root=ROOT):
    """Read-only exact pin check, shared by single-job and all-job entrypoints."""
    root, plan_path, freeze_path = Path(root), Path(plan_path), Path(freeze_path)
    freeze = json.loads(freeze_path.read_text())
    if freeze.get('status') != 'frozen_before_first_HADDOCK_cohort_pose' or not freeze.get('execution_authorized'):
        raise PairedCohortError('Root execution freeze is required')
    _checked(plan_path, freeze['prepared_plan']['sha256'])
    plan = json.loads(plan_path.read_text())
    if plan['ncores'] != 1 or plan['seed_cell_denominator'] != 108:
        raise PairedCohortError('Frozen CPU/denominator changed')
    for record in [plan['driver'], plan['proposal'], plan['assembly_protocol_policy'], plan['source_manifest'], plan['preflight_receipt'], plan['release_lock'], plan['CNS_binary'], plan['release_distribution'], plan['native_environment_freeze']] + plan['runtime_sources']:
        _checked(root / record['path'], record['sha256'])
    return plan


def execute_job(plan_path, freeze_path, job_id, root=ROOT):
    """One foreground job, one CPU, only after exact external root freeze."""
    root, plan_path, freeze_path = Path(root), Path(plan_path), Path(freeze_path)
    plan = verify_execution_freeze(plan_path, freeze_path, root)
    job = next(j for j in plan['jobs'] if j['job_id'] == job_id)
    if job['status'] != 'ready':
        return dict(job_id=job_id, status=job['status'], terminal=True, eligible_workflow_pool=True,
                    physical_models=0, no_native_execution=True)
    for record in [job['config']] + job['input_files']:
        _checked(root / record['path'], record['sha256'])
    directory = root / job['directory']; result_path = directory / 'execution.json'
    if result_path.exists() or (directory / 'run').exists():
        raise PairedCohortError('Refuse repeated or partial existing native run')
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1', MKL_NUM_THREADS='1')
    started = datetime.datetime.now(datetime.timezone.utc).isoformat(); start = time.monotonic()
    command = [str(root / 'tools/haddock3-2026.10.1/env/bin/haddock3'), 'run.cfg']
    with (directory / 'stdout.log').open('w') as stdout, (directory / 'stderr.log').open('w') as stderr:
        process = subprocess.run(command, cwd=directory, env=environment, stdout=stdout, stderr=stderr, check=False)
    stages = model_inventory(directory / 'run', job_id)
    initial = next((s for s in stages if s['stage'].endswith('_rigidbody')), None)
    final = next((s for s in stages if s['stage'].endswith('_flexref')), None)
    final_ranked = bool(final and final['actual_physical_model_files'] and
                        len(final['individual_score_order']) == final['actual_physical_model_files'])
    result = dict(schema_version=1, job_id=job_id, eligible_workflow_pool=True, terminal=True,
                  status='completed' if process.returncode == 0 and final_ranked else 'native_execution_failed',
                  process_return_code=process.returncode, started_utc=started, elapsed_seconds=time.monotonic() - start,
                  ncores=1, command=command, config=job['config'], prepared_plan=pin(plan_path, root),
                  root_freeze=pin(freeze_path, root), helper_source=plan['driver'], CNS_binary=plan['CNS_binary'],
                  stages=stages, no_accuracy_evaluation=True,
                  requested_initial_model_slots=20,
                  actual_initial_models=initial['actual_physical_model_files'] if initial else 0,
                  actual_refined_models=final['actual_physical_model_files'] if final else 0,
                  unfilled_initial_slots=20 - initial['actual_physical_model_files'] if initial else 20,
                  final_score_metadata_complete=final_ranked,
                  failed_model_slots_not_replaced=True, clustering_performed=False)
    _write(result_path, result)
    return result


def execute_all(plan_path, freeze_path, root=ROOT):
    """Sequential foreground traversal of the complete immutable denominator."""
    root, plan_path = Path(root), Path(plan_path)
    plan = verify_execution_freeze(plan_path, freeze_path, root)
    records = []
    for job in plan['jobs']:
        directory = root / job['directory'] if job.get('directory') else None
        existing = directory / 'execution.json' if directory else None
        if existing and existing.is_file():
            result = json.loads(existing.read_text())
            if not result.get('terminal') or result['prepared_plan']['sha256'] != hashlib.sha256(plan_path.read_bytes()).hexdigest():
                raise PairedCohortError('Existing execution does not match this frozen plan')
        else:
            result = execute_job(plan_path, freeze_path, job['job_id'], root)
        records.append(result)
        _write(plan_path.parent / 'execution_progress.json',
               dict(schema_version=1, prepared_plan=pin(plan_path, root), root_freeze=pin(Path(freeze_path), root),
                    seed_cell_denominator=108, terminal_cells=len(records), all_terminal=len(records) == 108,
                    no_accuracy_evaluation=True, jobs=records))
        print(job['job_id'], result['status'], flush=True)
    return records


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--plan', type=Path)
    parser.add_argument('--freeze', type=Path)
    parser.add_argument('--job-id')
    parser.add_argument('--run-all-foreground', action='store_true')
    args = parser.parse_args(argv)
    if args.prepare_only:
        plan = prepare_plan(args.root, args.output)
        print(json.dumps(dict(conditions=len(plan['conditions']), seed_cells=len(plan['jobs']),
                              ready=sum(j['status'] == 'ready' for j in plan['jobs'])), indent=2))
        return 0
    if args.run_all_foreground:
        if not args.plan or not args.freeze:
            parser.error('--run-all-foreground requires --plan and --freeze')
        records = execute_all(args.plan, args.freeze, args.root)
        return 0 if len(records) == 108 and all(r.get('terminal') for r in records) else 1
    if not all([args.plan, args.freeze, args.job_id]):
        parser.error('Use --prepare-only, --run-all-foreground, or all of --plan --freeze --job-id')
    result = execute_job(args.plan, args.freeze, args.job_id, args.root)
    print(json.dumps(dict(job_id=result['job_id'], status=result['status']), indent=2))
    return 0 if result['status'] == 'completed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
