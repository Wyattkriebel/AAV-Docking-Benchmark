"""Preserving orchestration recovery for the frozen HADDOCK main cohort.

The frozen method helper/plan/configs are never modified. A completed orphan is
adopted only from pinned raw evidence and with its lost return code left null.
Future calls normalize paths, use isolated process workers for resource records,
and require a separate coordinator recovery freeze before any native launch.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, FIRST_COMPLETED, wait
import datetime as dt
import gzip
import hashlib
import json
import math
import multiprocessing
import os
from pathlib import Path
import re
import resource
import time
import tomllib
from zoneinfo import ZoneInfo

import numpy as np

from . import haddock_paired_cohort as frozen
from .haddock_cohort_preparation import _checked, pin
from .haddock_method import pdb_atoms, read_cns_topology

ORPHAN_JOB = 'strict/2RDK/fpocket/box01/s104729'
PLAN_SHA = 'c447d1216d63cf5d19fb5a214a312cb35f2a67599f9c732a1a33276b3ca2c6f2'
HELPER_SHA = '2286023e47908c168b67a119ded26777919697f6ee682985b6938d0c9a7b5bf6'
ORIGINAL_FREEZE_SHA = 'f5c9a8503a6ebcb08cb9b5b6c3b1f698b36868c610a3816352829bf416cda904'
MAX_PARALLEL_JOBS = 6
OLD_PROGRESS = 'results/blind_pose/execution_progress.json'
RECOVERY_POLICY = 'benchmarks/haddock_orchestration_recovery_policy.json'


class RecoveryError(ValueError):
    pass


def normalize_paths(root, plan, original_freeze):
    root = Path(root).resolve()
    def resolve(path):
        path = Path(path)
        return (path if path.is_absolute() else root / path).resolve()
    return root, resolve(plan), resolve(original_freeze)


def _write_once(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        stream.write(json.dumps(data, indent=2) + '\n')


def _replace_json(path, data):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(data, indent=2) + '\n')
    temporary.replace(path)


def _text(path):
    path = Path(path)
    return gzip.open(path, 'rt').read() if path.suffix == '.gz' else path.read_text()


def _native_context(root, plan, original_freeze):
    root, plan, original_freeze = normalize_paths(root, plan, original_freeze)
    _checked(plan, PLAN_SHA)
    _checked(original_freeze, ORIGINAL_FREEZE_SHA)
    _checked(root / 'aav_benchmark/haddock_paired_cohort.py', HELPER_SHA)
    specification = frozen.verify_execution_freeze(plan, original_freeze, root)
    return root, plan, original_freeze, specification


def _atom_key(atom):
    return tuple(atom[k] for k in ['chain', 'resnum', 'insertion', 'resname', 'atom_name', 'alt', 'element'])


def coordinate_identity(before, after):
    """Allow header/TER normalization, require every physical atom unchanged."""
    left = {_atom_key(a): a['xyz_A'] for a in pdb_atoms(before)}
    right = {_atom_key(a): a['xyz_A'] for a in pdb_atoms(after)}
    if set(left) != set(right):
        raise RecoveryError('Native sanitized input changed atom identity')
    maximum = max((float(np.linalg.norm(np.array(left[k]) - right[k])) for k in left), default=0)
    if maximum != 0:
        raise RecoveryError('Native sanitized input changed coordinates')
    return dict(atom_count=len(left), exact_atom_identity=True, maximum_coordinate_change_A=maximum,
                byte_identical=Path(before).read_bytes() == Path(after).read_bytes())


def completion_log_proof(stdout, original_freeze, timezone='America/Los_Angeles'):
    text = Path(stdout).read_text()
    durations = re.findall(r'This HADDOCK3 run took:\s*(\d+) minutes? and (\d+) seconds', text)
    if len(durations) != 1 or len(re.findall(r'\bcli INFO\] Finished at ', text)) != 1:
        raise RecoveryError('Missing or ambiguous normal native completion/runtime marker')
    if 'Traceback (most recent call last)' in text:
        raise RecoveryError('Native stdout contains a traceback')
    matches = re.findall(r'^\[(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d,\d+)\s', text, flags=re.M)
    if not matches:
        raise RecoveryError('Missing timestamped native start message')
    first = dt.datetime.strptime(matches[0], '%Y-%m-%d %H:%M:%S,%f').replace(tzinfo=ZoneInfo(timezone))
    original = json.loads(Path(original_freeze).read_text())
    freeze_time = dt.datetime.fromisoformat(original['frozen_at_UTC'])
    if first.astimezone(dt.timezone.utc) < freeze_time:
        raise RecoveryError('Timestamped native message predates original freeze')
    minutes, seconds = map(int, durations[0])
    return dict(normal_finished_marker=True, elapsed_seconds=minutes * 60 + seconds,
                elapsed_source='explicit integer-second HADDOCK native stdout runtime report',
                measured_external_elapsed_seconds=None, process_CPU_seconds=None,
                first_native_log_timestamp=matches[0], native_log_timezone_context=timezone,
                first_native_log_UTC_derived=first.astimezone(dt.timezone.utc).isoformat(),
                process_return_code=None, return_code_reason='Not captured before original post-native Python receipt serialization exception; completion markers do not establish exit0.')


def completed_output_proof(directory, job_id, job):
    """Require all20 initial/refined products, exact lineage and finite records."""
    directory = Path(directory)
    stages = frozen.model_inventory(directory / 'run', job_id)
    body = next((s for s in stages if s['stage'].endswith('_rigidbody')), None)
    refined = next((s for s in stages if s['stage'].endswith('_flexref')), None)
    if not body or not refined:
        raise RecoveryError('Missing complete native stages')
    topology_dir = next((p for p in (directory / 'run').glob('*_topoaa') if p.is_dir()), None)
    if topology_dir is None:
        raise RecoveryError('Missing actual topology stage')
    expected = set()
    psf_counts = {}
    for molecule in ['receptor', 'glycan']:
        pdb = topology_dir / (molecule + '_haddock.pdb')
        psf = topology_dir / (molecule + '_haddock.psf')
        atoms = pdb_atoms(pdb)
        topology = read_cns_topology(psf)
        if len(atoms) != topology['atom_count']:
            raise RecoveryError('Prepared coordinates and topology atom counts differ')
        expected.update(_atom_key(a) for a in atoms)
        psf_counts[molecule] = topology['atom_count']
    lineage = []
    body_io = json.loads((directory / 'run' / body['stage'] / 'io.json').read_text())
    refined_io = json.loads((directory / 'run' / refined['stage'] / 'io.json').read_text())
    body_metadata = {r['file_name']: r for r in body_io['output']}
    refined_metadata = {r['file_name']: r for r in refined_io['output']}
    body_names = set(body_metadata)
    input_names = {r['file_name'] for r in refined_io['input']}
    if input_names != body_names:
        raise RecoveryError('Refinement input pool differs from all generated initial models')
    for stage in [body, refined]:
        if stage['actual_physical_model_files'] != 20 or len(stage['individual_score_order']) != 20:
            raise RecoveryError('Adoption requires complete20 physical and score-ranked models per stage')
        if {m['original_model_number'] for m in stage['models']} != set(range(1, 21)):
            raise RecoveryError('Missing or duplicate native model slots')
        if not (directory / 'run' / stage['stage'] / 'params.cfg').is_file():
            raise RecoveryError('Missing effective native stage parameters')
        for model in stage['models']:
            path = directory / 'run' / stage['stage'] / model['file_name']
            atoms = pdb_atoms(path)
            if {_atom_key(a) for a in atoms} != expected or len(atoms) != len(expected):
                raise RecoveryError('Model atom inventory differs from actual topology')
            if not all(np.isfinite(a['xyz_A']).all() for a in atoms):
                raise RecoveryError('Nonfinite native model coordinates')
            if not model['energies'] or not all(isinstance(v, (int, float)) and math.isfinite(v) for v in model['energies'].values()):
                raise RecoveryError('Missing/nonfinite native method energy inventory')
    if {r['seed'] for r in body_metadata.values()} != set(range(job['seed'] + 1, job['seed'] + 21)):
        raise RecoveryError('Actual initial model seeds differ from frozen seed range')
    for name, row in refined_metadata.items():
        parent = row.get('ori_name')
        if parent not in body_metadata or row.get('seed') != body_metadata[parent]['seed']:
            raise RecoveryError('Missing or changed refinement lineage/seed')
        lineage.append(dict(refined_model=name, initial_model=parent, actual_CNS_seed=row['seed']))
    if len({r['initial_model'] for r in lineage}) != 20:
        raise RecoveryError('Refined pool duplicates or omits an initial parent')
    logs = []
    for stage in [body, refined]:
        paths = sorted((directory / 'run' / stage['stage']).glob('*.out*'))
        if len(paths) != 20:
            raise RecoveryError('Missing actual CNS model logs')
        for path in paths:
            text = _text(path)
            if 'Program stopped at:' not in text or 'CNSsolve>stop' not in text or re.search(r'%[A-Za-z0-9_-]*ERR|^\s*(?:FATAL|ABORT)', text, re.M | re.I):
                raise RecoveryError('CNS model log lacks normal stop or contains fatal markers')
            logs.append(path)
    return dict(stages=stages, actual_initial_models=20, actual_refined_models=20,
                topology_atom_counts=psf_counts, full_model_atom_count=len(expected),
                complete_model_score_inventory=True, complete_finite_energy_inventory=True,
                refinement_lineage=sorted(lineage, key=lambda r: r['refined_model']), normal_CNS_model_logs=len(logs))


def audit_orphan(root, plan_path, original_freeze, job_id=ORPHAN_JOB):
    root, plan_path, original_freeze, plan = _native_context(root, plan_path, original_freeze)
    if job_id != ORPHAN_JOB:
        raise RecoveryError('Only the explicitly reviewed first orphan is eligible for this adoption')
    job = next(j for j in plan['jobs'] if j['job_id'] == job_id)
    directory = root / job['directory']
    if (directory / 'execution.json').exists():
        raise RecoveryError('Orphan already has an execution receipt; do not overwrite')
    for record in [job['config']] + job['input_files']:
        _checked(root / record['path'], record['sha256'])
    raw_config = directory / 'run/data/configurations/raw_input.toml'
    if raw_config.read_bytes() != (root / job['config']['path']).read_bytes():
        raise RecoveryError('Raw native config differs from frozen config bytes')
    config = tomllib.loads(raw_config.read_text())
    identities = {}
    for molecule in ['receptor', 'glycan']:
        identities[molecule] = coordinate_identity(directory.parent / 'data' / (molecule + '.pdb'),
                                                   directory / 'run/data/0_topoaa' / (molecule + '.pdb'))
    air = next(r for r in job['input_files'] if Path(r['path']).name == 'box_air.tbl')
    for stage in ['1_rigidbody', '2_flexref']:
        _checked(directory / 'run/data' / stage / 'box_air.tbl', air['sha256'])
    for stage in ['topoaa', 'rigidbody', 'flexref']:
        stage_dir = next(p for p in (directory / 'run').glob('*_' + stage) if p.is_dir())
        actual = tomllib.loads((stage_dir / 'params.cfg').read_text())[stage]
        for key, value in config[stage].items():
            if key == 'ambig_fname':
                continue
            if actual.get(key) != value:
                raise RecoveryError('Effective native parameter differs from frozen explicit config: ' + stage + '.' + key)
    log = completion_log_proof(directory / 'stdout.log', original_freeze)
    outputs = completed_output_proof(directory, job_id, job)
    failure_log = plan_path.parent / 'execution_driver.log'
    failure_text = failure_log.read_text()
    if 'ValueError:' not in failure_text or 'path.relative_to(root)' not in failure_text:
        raise RecoveryError('Original post-native receipt serialization exception not retained')
    raw_files = [pin(p, root) for p in sorted(directory.rglob('*')) if p.is_file()]
    proof = dict(schema_version=1, job_id=job_id, source_only_recovery=True, no_native_rerun=True,
                 no_accuracy_evaluation=True, frozen_plan=pin(plan_path, root), frozen_helper=plan['driver'],
                 original_execution_freeze=pin(original_freeze, root), frozen_config=job['config'],
                 original_failed_driver_log=pin(failure_log, root), input_preprocessing_identity=identities,
                 native_completion=log, output_proof=outputs, raw_artifacts=raw_files,
                 posthoc_hash_timing='Output hashes are anchored now after the metadata failure; they are not claimed to have been recorded during the lost subprocess telemetry.',
                 process_return_code=None, process_CPU_seconds=None, external_measured_wall_seconds=None)
    return proof


def _existing_receipts(root, plan):
    rows = []
    for job in plan['jobs']:
        if not job.get('directory'):
            continue
        path = root / job['directory'] / 'execution.json'
        if path.exists():
            row = json.loads(path.read_text())
            if not row.get('terminal') or row.get('job_id') != job['job_id'] or row.get('prepared_plan', {}).get('sha256') != PLAN_SHA:
                raise RecoveryError('Existing terminal receipt does not match original frozen job/plan')
            if row.get('config') != job['config']:
                raise RecoveryError('Existing terminal receipt names a different config')
            rows.append(pin(path, root))
    return rows


def prepare_recovery(root, plan_path, original_freeze, output):
    root, plan_path, original_freeze, plan = _native_context(root, plan_path, original_freeze)
    output = (root / output if not Path(output).is_absolute() else Path(output)).resolve()
    proof = audit_orphan(root, plan_path, original_freeze)
    output.mkdir(parents=True, exist_ok=True)
    _write_once(output / 'orphan_completion_proof.json', proof)
    recovery = dict(schema_version=1, status='prepared_waiting_for_root_recovery_freeze',
                    original_plan=pin(plan_path, root), original_execution_freeze=pin(original_freeze, root),
                    frozen_helper=plan['driver'], recovery_wrapper=pin(Path(__file__), root),
                    recovery_policy=pin(root / RECOVERY_POLICY, root),
                    orphan_completion_proof=pin(output / 'orphan_completion_proof.json', root),
                    orphan_job_id=ORPHAN_JOB, existing_terminal_receipts=_existing_receipts(root, plan),
                    seed_cell_denominator=108, original_ready_jobs=39, recovered_completed_jobs=1,
                    remaining_unique_native_jobs=38, native_cpus_per_job=1,
                    allowed_parallel_jobs_maximum=MAX_PARALLEL_JOBS,
                    native_model_parameters_changed=False, no_native_rerun=True, no_accuracy_evaluation=True)
    _write_once(output / 'recovery_plan.json', recovery)
    return recovery


def verify_recovery_authorization(root, recovery_plan_path, recovery_freeze_path):
    root = Path(root).resolve()
    recovery_plan_path = (root / recovery_plan_path if not Path(recovery_plan_path).is_absolute() else Path(recovery_plan_path)).resolve()
    recovery_freeze_path = (root / recovery_freeze_path if not Path(recovery_freeze_path).is_absolute() else Path(recovery_freeze_path)).resolve()
    authorization = json.loads(recovery_freeze_path.read_text())
    if authorization.get('status') != 'frozen_before_HADDOCK_orchestration_resumption' or not authorization.get('execution_authorized'):
        raise RecoveryError('Separate coordinator recovery freeze required')
    _checked(recovery_plan_path, authorization['recovery_plan']['sha256'])
    recovery = json.loads(recovery_plan_path.read_text())
    _checked(root / recovery['recovery_wrapper']['path'], recovery['recovery_wrapper']['sha256'])
    _checked(root / recovery['recovery_policy']['path'], recovery['recovery_policy']['sha256'])
    _checked(root / recovery['orphan_completion_proof']['path'], recovery['orphan_completion_proof']['sha256'])
    proof = json.loads((root / recovery['orphan_completion_proof']['path']).read_text())
    for record in proof['raw_artifacts']:
        _checked(root / record['path'], record['sha256'])
    cap = authorization['maximum_parallel_jobs']
    if not isinstance(cap, int) or isinstance(cap, bool) or not 1 <= cap <= MAX_PARALLEL_JOBS:
        raise RecoveryError('Invalid authorized parallel-job cap')
    if authorization.get('native_cpus_per_job') != 1:
        raise RecoveryError('Per-job native CPU remains one')
    for record in recovery['existing_terminal_receipts']:
        _checked(root / record['path'], record['sha256'])
    root, plan_path, original_freeze, plan = _native_context(root, recovery['original_plan']['path'], recovery['original_execution_freeze']['path'])
    return root, recovery_plan_path, recovery_freeze_path, authorization, recovery, plan_path, original_freeze, plan


def adopt_completed_outputs(root, recovery_plan_path, recovery_freeze_path):
    root, recovery_plan_path, recovery_freeze_path, authorization, recovery, plan_path, original_freeze, plan = verify_recovery_authorization(root, recovery_plan_path, recovery_freeze_path)
    proof_path = root / recovery['orphan_completion_proof']['path']
    proof = json.loads(proof_path.read_text())
    job = next(j for j in plan['jobs'] if j['job_id'] == ORPHAN_JOB)
    directory = root / job['directory']; execution = directory / 'execution.json'
    if execution.exists():
        existing = json.loads(execution.read_text())
        if existing.get('status') != 'complete_outputs_recovered' or existing.get('process_return_code') is not None:
            raise RecoveryError('Refuse to overwrite or reinterpret existing receipt')
        _checked(root / existing['recovery_receipt']['path'], existing['recovery_receipt']['sha256'])
        return existing
    actual = audit_orphan(root, plan_path, original_freeze)
    if actual != proof:
        raise RecoveryError('Orphan evidence changed after recovery audit')
    receipt_path = recovery_plan_path.parent / 'adoption_receipt.json'
    receipt = dict(schema_version=1, status='verified_completed_outputs_adopted_without_rerun', job_id=ORPHAN_JOB,
                   completion_proof=pin(proof_path, root), recovery_plan=pin(recovery_plan_path, root),
                   recovery_freeze=pin(recovery_freeze_path, root), original_plan=pin(plan_path, root),
                   original_driver_failure=proof['original_failed_driver_log'], native_rerun=False,
                   original_native_return_code=None, original_process_CPU_seconds=None,
                   native_self_reported_elapsed_seconds=proof['native_completion']['elapsed_seconds'],
                   no_accuracy_evaluation=True)
    if receipt_path.exists():
        if json.loads(receipt_path.read_text()) != receipt:
            raise RecoveryError('Existing adoption receipt changed')
    else:
        _write_once(receipt_path, receipt)
    result = dict(schema_version=1, job_id=ORPHAN_JOB, eligible_workflow_pool=True, terminal=True,
                  status='complete_outputs_recovered', process_return_code=None,
                  process_return_code_provenance=proof['native_completion']['return_code_reason'],
                  started_utc=None, elapsed_seconds=proof['native_completion']['elapsed_seconds'],
                  elapsed_seconds_provenance=proof['native_completion']['elapsed_source'],
                  external_measured_wall_seconds=None, process_CPU_seconds=None, ncores=1,
                  config=job['config'], prepared_plan=pin(plan_path, root), root_freeze=pin(original_freeze, root),
                  execution_recovery_freeze=pin(recovery_freeze_path, root), helper_source=plan['driver'],
                  recovery_wrapper_source=recovery['recovery_wrapper'], recovery_receipt=pin(receipt_path, root),
                  CNS_binary=plan['CNS_binary'], stages=proof['output_proof']['stages'],
                  requested_initial_model_slots=20, actual_initial_models=20, actual_refined_models=20,
                  unfilled_initial_slots=0, final_score_metadata_complete=True,
                  failed_model_slots_not_replaced=True, clustering_performed=False,
                  native_rerun=False, no_accuracy_evaluation=True)
    _write_once(execution, result)
    return result


def resource_snapshot():
    def values(kind):
        usage = resource.getrusage(kind)
        return dict(user_seconds=usage.ru_utime, system_seconds=usage.ru_stime,
                    maximum_rss=usage.ru_maxrss)
    return dict(worker=values(resource.RUSAGE_SELF), waited_children=values(resource.RUSAGE_CHILDREN))


def measured_usage(before, after, wall):
    components = {}
    for scope in ['worker', 'waited_children']:
        components[scope] = {key: after[scope][key] - before[scope][key] for key in ['user_seconds', 'system_seconds']}
        if min(components[scope].values()) < 0:
            raise RecoveryError('OS resource counter went backwards')
    return dict(wall_seconds=wall, CPU_seconds=sum(sum(c.values()) for c in components.values()),
                components=components, before=before, after=after,
                scope='Single isolated job worker call plus its waited-on HADDOCK subprocess descendants; worker startup before the call is excluded. RSS is reported as raw OS counters, not subtracted or relabeled as per-job peak.',
                measurement='time.perf_counter and resource.getrusage(RUSAGE_SELF/RUSAGE_CHILDREN)')


def measured_native_worker(root, plan_path, original_freeze, job_id):
    """Exactly one original native job in an isolated spawned process."""
    root, plan_path, original_freeze = normalize_paths(root, plan_path, original_freeze)
    plan = json.loads(plan_path.read_text()); job = next(j for j in plan['jobs'] if j['job_id'] == job_id)
    directory = root / job['directory']; before = resource_snapshot(); start = time.perf_counter(); failure = None
    try:
        result = frozen.execute_job(plan_path, original_freeze, job_id, root)
    except Exception as error:
        failure = dict(exception_type=type(error).__name__, reason=str(error))
        try:
            partial = frozen.model_inventory(directory / 'run', job_id) if (directory / 'run').is_dir() else []
        except Exception as inventory_error:
            partial = []
            failure['inventory_error'] = dict(exception_type=type(inventory_error).__name__, reason=str(inventory_error))
        result = dict(schema_version=1, job_id=job_id, status='orchestration_failed', terminal=True,
                      eligible_workflow_pool=True, process_return_code=None, no_accuracy_evaluation=True,
                      prepared_plan=pin(plan_path, root), config=job['config'],
                      retained_partial_stages=partial,
                      failure=failure, native_retry=False)
        _write_once(directory / 'wrapper_failure.json', result)
    finally:
        after = resource_snapshot(); wall = time.perf_counter() - start
        sidecar = dict(schema_version=1, job_id=job_id, native_cpus_per_job=1,
                       worker_pid=os.getpid(), measured_usage=measured_usage(before, after, wall),
                       original_plan=pin(plan_path, root), frozen_helper=plan['driver'],
                       recovery_wrapper=pin(Path(__file__), root), failure=failure)
        _write_once(directory / 'measured_resources.json', sidecar)
    return dict(result=result, resource_usage=pin(directory / 'measured_resources.json', root),
                execution_receipt=pin(directory / ('wrapper_failure.json' if failure else 'execution.json'), root))


def read_control(path, cap):
    control = json.loads(Path(path).read_text())
    if set(control) != {'active_job_limit', 'stop_after_running'}:
        raise RecoveryError('Concurrency control has unexpected keys')
    number = control['active_job_limit']
    if not isinstance(number, int) or isinstance(number, bool) or not 1 <= number <= cap or not isinstance(control['stop_after_running'], bool):
        raise RecoveryError('Concurrency control exceeds authorized bounds')
    return control


def workflow_concurrency_limit(progress, maximum):
    """Predeclared resource rule; inspect workflow completion, never accuracy."""
    complete = (progress.get('scheduled_seed_runs') == 156 and
                progress.get('canonical_terminal_seed_runs') == 156 and
                progress.get('running_seed_runs') == 0 and progress.get('status') == 'complete')
    return min(maximum, 6 if complete else 2)


def _run_resume_pool_impl(root, recovery_plan_path, recovery_freeze_path, control_path):
    context = verify_recovery_authorization(root, recovery_plan_path, recovery_freeze_path)
    root, recovery_plan_path, recovery_freeze_path, authorization, recovery, plan_path, original_freeze, plan = context
    control_path = (root / control_path if not Path(control_path).is_absolute() else Path(control_path)).resolve()
    if str(control_path.relative_to(root)) != authorization['concurrency_control_path']:
        raise RecoveryError('Unexpected concurrency control path')
    cap = authorization['maximum_parallel_jobs']; read_control(control_path, cap)
    required_rule = dict(progress_path=OLD_PROGRESS, scheduled_seed_runs=156,
                         before_complete_maximum=2, after_complete_maximum=6,
                         required_final_status='complete', only_workflow_status_used=True)
    if authorization.get('operational_concurrency_rule') != required_rule:
        raise RecoveryError('Missing exact predeclared old-workflow resource rule')
    adopt_completed_outputs(root, recovery_plan_path, recovery_freeze_path)
    terminal, pending, immutable = {}, [], {}
    for job in plan['jobs']:
        if job['status'] != 'ready':
            terminal[job['job_id']] = dict(result=dict(job_id=job['job_id'], status=job['status'], terminal=True,
                eligible_workflow_pool=True, physical_models=0, no_native_execution=True))
            continue
        directory = root / job['directory']; existing = directory / 'execution.json'; failure = directory / 'wrapper_failure.json'
        path = existing if existing.exists() else failure if failure.exists() else None
        if path:
            row = json.loads(path.read_text())
            if not row.get('terminal') or row.get('job_id') != job['job_id'] or row.get('prepared_plan', {}).get('sha256') != PLAN_SHA:
                raise RecoveryError('Existing receipt not tied to frozen job')
            terminal[job['job_id']] = dict(result=row, execution_receipt=pin(path, root))
            if (directory / 'measured_resources.json').exists():
                terminal[job['job_id']]['resource_usage'] = pin(directory / 'measured_resources.json', root)
            immutable[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        elif (directory / 'run').exists():
            raise RecoveryError('Unreviewed partial native run: do not rerun')
        else:
            pending.append(job['job_id'])
    progress_path = plan_path.parent / 'execution_progress.json'
    history = []
    if progress_path.exists():
        previous = json.loads(progress_path.read_text())
        if previous.get('prepared_plan', {}).get('sha256') != PLAN_SHA:
            raise RecoveryError('Existing progress belongs to a different frozen plan')
        history = previous.get('orchestration_control_history', [])
    last_control = None
    def progress():
        for name, sha in immutable.items():
            _checked(Path(name), sha)
        ordered = []
        for job in plan['jobs']:
            if job['job_id'] in terminal:
                wrapped = terminal[job['job_id']]; row = dict(wrapped['result'])
                for key in ['execution_receipt', 'resource_usage']:
                    if key in wrapped:
                        row[key] = wrapped[key]
                ordered.append(row)
        _replace_json(progress_path, dict(schema_version=2, prepared_plan=pin(plan_path, root), root_freeze=pin(original_freeze, root),
            execution_recovery_freeze=pin(recovery_freeze_path, root), recovery_plan=pin(recovery_plan_path, root),
            seed_cell_denominator=108, terminal_cells=len(ordered), all_terminal=len(ordered) == 108,
            original_native_cpus_per_job=1, authorized_parallel_jobs_maximum=cap,
            orchestration_control_history=history, no_accuracy_evaluation=True, jobs=ordered))
    progress()
    with ProcessPoolExecutor(max_workers=cap, mp_context=multiprocessing.get_context('spawn'), max_tasks_per_child=1) as executor:
        active = {}
        while pending or active:
            control = read_control(control_path, cap)
            try:
                resource_bytes = (root / OLD_PROGRESS).read_bytes()
                resource_status = json.loads(resource_bytes)
                status_read_valid = True
            except (FileNotFoundError, json.JSONDecodeError):
                resource_bytes = b'{}\n'
                resource_status = {}
                status_read_valid = False
            operational_limit = workflow_concurrency_limit(resource_status, cap)
            effective_limit = min(control['active_job_limit'], operational_limit)
            effective_control = dict(control, effective_active_job_limit=effective_limit)
            if effective_control != last_control:
                snapshots = recovery_plan_path.parent / 'resource_status_snapshots'
                snapshots.mkdir(exist_ok=True)
                index = max([int(p.stem) for p in snapshots.glob('*.json') if p.stem.isdigit()] + [0]) + 1
                snapshot = snapshots / f'{index:04d}.json'
                with snapshot.open('xb') as stream:
                    stream.write(resource_bytes)
                history.append(dict(observed_utc=dt.datetime.now(dt.timezone.utc).isoformat(), control=control,
                                    effective_active_job_limit=effective_limit,
                                    control_source=pin(control_path, root),
                                    historical_workflow_status_snapshot=pin(snapshot, root),
                                    status_read_valid=status_read_valid,
                                    missing_or_partial_status_means_fail_closed_to_two=not status_read_valid,
                                    original_status_path=OLD_PROGRESS, no_accuracy_fields_consulted=True))
                last_control = effective_control; progress()
            while pending and not control['stop_after_running'] and len(active) < effective_limit:
                job_id = pending.pop(0)
                future = executor.submit(measured_native_worker, str(root), str(plan_path), str(original_freeze), job_id)
                active[future] = job_id
                print('START', job_id, 'one_native_CPU', flush=True)
            if not active:
                break
            done, _ = wait(active, timeout=1, return_when=FIRST_COMPLETED)
            for future in done:
                job_id = active.pop(future)
                try:
                    wrapped = future.result()
                except Exception as error:
                    raise RecoveryError('Isolated worker failed without a verified terminal receipt; no rerun: ' + job_id) from error
                terminal[job_id] = wrapped
                path = root / wrapped['execution_receipt']['path']
                immutable[str(path)] = wrapped['execution_receipt']['sha256']
                progress(); print('TERMINAL', job_id, wrapped['result']['status'], flush=True)
    progress()
    return dict(all_terminal=len(terminal) == 108, terminal_cells=len(terminal), pending_cells=len(pending))


def run_resume_pool(root, recovery_plan_path, recovery_freeze_path, control_path):
    context = verify_recovery_authorization(root, recovery_plan_path, recovery_freeze_path)
    root, recovery_plan_path = context[:2]
    lock = recovery_plan_path.parent / 'coordinator.lock'
    try:
        with lock.open('x') as stream:
            stream.write(json.dumps(dict(pid=os.getpid(), recovery_plan=pin(recovery_plan_path, root),
                                        created_utc=dt.datetime.now(dt.timezone.utc).isoformat())) + '\n')
    except FileExistsError as error:
        raise RecoveryError('Another or unreviewed prior coordinator lock exists; do not submit duplicate jobs') from error
    try:
        return _run_resume_pool_impl(root, recovery_plan_path, recovery_freeze_path, control_path)
    finally:
        if lock.exists() and json.loads(lock.read_text()).get('pid') == os.getpid():
            lock.unlink()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=frozen.ROOT)
    parser.add_argument('--plan', type=Path, default=frozen.DEFAULT_OUTPUT / 'prepared_plan.json')
    parser.add_argument('--original-freeze', type=Path, default=frozen.DEFAULT_OUTPUT / 'execution_freeze.json')
    parser.add_argument('--prepare-recovery', action='store_true')
    parser.add_argument('--recovery-directory', type=Path, default=frozen.DEFAULT_OUTPUT / 'orchestration_recovery')
    parser.add_argument('--recovery-plan', type=Path)
    parser.add_argument('--recovery-freeze', type=Path)
    parser.add_argument('--control', type=Path)
    args = parser.parse_args(argv)
    if args.prepare_recovery:
        record = prepare_recovery(args.root, args.plan, args.original_freeze, args.recovery_directory)
        print(json.dumps(dict(status=record['status'], remaining_unique_native_jobs=38), indent=2)); return 0
    if not all([args.recovery_plan, args.recovery_freeze, args.control]):
        parser.error('Use --prepare-recovery or --recovery-plan --recovery-freeze --control')
    result = run_resume_pool(args.root, args.recovery_plan, args.recovery_freeze, args.control)
    print(json.dumps(result, indent=2)); return 0


if __name__ == '__main__':
    raise SystemExit(main())
