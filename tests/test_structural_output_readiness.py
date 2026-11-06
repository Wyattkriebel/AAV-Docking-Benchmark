"""Complete synthetic TEMP metadata cohorts; never inspect actual workflows."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import socket
import shutil
import subprocess
import sys
from unittest.mock import patch

import pytest

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'results/structural_pose_output_foundation/validate_readiness.py'
SPEC = importlib.util.spec_from_file_location('synthetic_structural_readiness', SOURCE)
validator = importlib.util.module_from_spec(SPEC)
previous_bytecode = sys.dont_write_bytecode
sys.dont_write_bytecode = True
SPEC.loader.exec_module(validator)
HELPER_SPEC = importlib.util.spec_from_file_location('readiness_packaging_test_helpers', ROOT / 'tests/test_structural_output_packager.py')
helpers = importlib.util.module_from_spec(HELPER_SPEC)
HELPER_SPEC.loader.exec_module(helpers)
sys.dont_write_bytecode = previous_bytecode
pkg = helpers.pkg


@pytest.fixture
def tmp_path(tmp_path_factory):
    return tmp_path_factory.mktemp('structural-readiness').resolve()


@pytest.fixture(autouse=True)
def no_bytecode(monkeypatch):
    monkeypatch.setattr(sys, 'dont_write_bytecode', True)


def put(root, name, raw):
    return helpers.put(root, name, raw)


def jput(root, name, value):
    return put(root, name, validator.canonical_bytes(value) + b'\n')


def pin(root, name):
    raw = (root / name).read_bytes()
    return {'path': name, 'bytes': len(raw), 'sha256': pkg.digest(raw)}


def model_blob(count=2):
    return b''.join(f'MODEL {rank}\n'.encode() + b'REMARK arbitrary energy text is opaque\nATOM coordinate fields are intentionally not numerical\nENDMDL\n'
                    for rank in range(1, count + 1))


def streams(**extra):
    value = {'raw_stdout': '{"opaque_energy":1.234567890123456789e+9999}\n',
             'raw_stderr': 'unknown diagnostic text\n', 'recognition_log': 'opaque recognition stream\n'}
    value.update({name + '_sha256': pkg.digest(value[name].encode()) for name in list(value)})
    return {**value, **extra}


def configuration():
    return {'seeds': validator.SEEDS, 'cpu': 2, 'exhaustiveness': 16, 'num_modes': 20,
            'energy_range_kcal_mol': 10, 'min_rmsd_A': validator.OpaqueFloat('0.5')}


def modern_condition(root, row, config, config_pin, index_sha, *, assembly=False, fail_one=False):
    receipt_path = str(Path(row['modern_runs']).parent.parent / 'execution.json')
    binding = 'prepared_conditions_sha256' if assembly else 'frozen_input_manifest_sha256'
    receipt = {'condition_id': row['condition_id'], 'configuration_sha256': config_pin['sha256'],
               'source_condition': row['source_condition'], binding: index_sha}
    outcomes = {}
    if row['status'] != 'prepared':
        receipt.update(status='prepared_input_failure', native_search_started=False, eligible_seed_slots=validator.SEEDS,
                       no_runs_or_poses_invented=True, failure_reason='synthetic source gate')
        jput(root, receipt_path, receipt)
        return {s: {'status': 'prepared_input_failure', 'record': receipt} for s in validator.SEEDS}, receipt_path, receipt
    prepared = row['modern_vina_prepared_directory']; prep = validator.load_json((root / row['preparation']['path']).read_bytes())
    directory = str(Path(row['modern_runs']).parent); runs = []
    jput(root, directory + '/configuration.json', config)
    for index, seed in enumerate(validator.SEEDS):
        status = 'failed' if fail_one and index == 0 else 'ok'
        pose_file, log_file = f'neighborhood_{seed}.pdbqt', f'neighborhood_{seed}.log'
        pose = put(root, directory + '/' + pose_file, model_blob() if status == 'ok' else b'opaque partial failed output\n')
        put(root, directory + '/' + log_file, b'opaque native output; no score interpretation\n')
        fingerprint = {'executable_sha256': 'e' * 64, 'config_sha256': config_pin['sha256'],
                       'preparation_sha256': row['preparation']['sha256'], 'ligand_sha256': row['ligand']['sha256'],
                       'receptor_sha256': row['receptor']['sha256'], 'receptor_pdb_sha256': pin(root, prepared + '/neighborhood.pdb')['sha256'],
                       'reference_sha256': pin(root, prepared + '/reference.json')['sha256'], 'full_assembly_sha256': None}
        command = ['/synthetic/vina', '--receptor', prepared + '/neighborhood.pdbqt', '--ligand', prepared + '/ligand.pdbqt',
                   '--out', directory + '/' + pose_file, '--scoring', 'vina', '--cpu', '2', '--seed', str(seed),
                   '--exhaustiveness', '16', '--num_modes', '20', '--min_rmsd', '0.5', '--energy_range', '10']
        for axis, center, size in zip('xyz', prep['box_center'], prep['box_size']):
            command += ['--center_' + axis, str(center), '--size_' + axis, str(size)]
        record = {'seed': seed, 'receptor': 'neighborhood', 'case_id': row['case_id'], 'status': status,
                  'information_regime': 'true-site-local' if row['method'] == 'oracle' else 'predicted-site-local',
                  'vina_version': 'AutoDock Vina v1.2.7', 'fingerprint': fingerprint, 'command': command,
                  'pose_file': pose_file, 'log_file': log_file, 'exit_code': 0 if status == 'ok' else 1,
                  'wall_seconds': 1, 'cpu_seconds': 1, 'child_peak_rss_bytes': 100,
                  'opaque_scientific_number': validator.OpaqueFloat('1.234567890123456789e+9999')}
        if status == 'ok': record.update(returned_poses=2, pose_sha256=pose['sha256'])
        else: record['error'] = 'synthetic failed native invocation'
        runs.append(record); outcomes[seed] = {'status': status, 'record': record, 'pose_path': pose['path'], 'models': 2 if status == 'ok' else 0}
    run_pin = jput(root, row['modern_runs'], runs)
    receipt.update(status='complete', native_search_started=True, terminal_seed_runs=3,
                   seed_statuses={str(r['seed']): r['status'] for r in runs}, runs_path=row['modern_runs'], runs_sha256=run_pin['sha256'])
    jput(root, receipt_path, receipt)
    return outcomes, receipt_path, receipt


def source_row(root, profile, case, arm, prepared, config_pin, *, assembly=False):
    base = validator.ASSEMBLY if assembly else validator.MAIN
    condition_id = f'{profile}/{case}/{arm}/box01'; directory = f'{base}/{profile}/{case}/{arm}/prepared'
    condition = {'case_id': case, 'method': arm, 'profile': profile, 'configuration_sha256': config_pin['sha256'],
                 'status': 'prepared' if prepared else 'preparation_failed', 'jobs': []}
    row = {'condition_id': condition_id, 'case_id': case, 'receptor_profile': profile, 'method': arm,
           'eligible_workflow_condition': True, 'status': 'prepared' if prepared else 'preparation_failed',
           'modern_runs': f'{base}/{profile}/{case}/{arm}/vina/box01/runs.json'}
    if prepared:
        box = directory + '/box01'
        row['ligand'] = put(root, box + '/ligand.pdbqt', b'opaque independent ligand\n')
        row['receptor'] = put(root, box + '/neighborhood.pdbqt', b'opaque prepared receptor\n')
        put(root, box + '/neighborhood.pdb', b'opaque receptor PDB\n')
        jput(root, box + '/reference.json', {'names': ['C1'], 'pose_reference_eligible': False})
        row['canonical_atom_map'] = jput(root, box + '/ligand_atom_map.json', {'atoms': [{'name': 'C1', 'source_key': ['B', 1, '', 'MAN', 'C1']}]})
        prep = {'case_id': case, 'profile': profile, 'method': arm, 'configuration_sha256': config_pin['sha256'],
                'ligand_file': 'ligand.pdbqt', 'box_center': [1, 2, 3], 'box_size': [24, 24, 24],
                'execution_information_regime': 'true-site-local' if arm == 'oracle' else 'predicted-site-local',
                'receptors': [{'id': 'neighborhood', 'state': 'holo', 'pdb_file': 'neighborhood.pdb', 'pdbqt_file': 'neighborhood.pdbqt'}]}
        row['preparation'] = jput(root, box + '/preparation.json', prep)
        row.update(candidate_class_rank=1, candidate_original_rank=1, modern_vina_prepared_directory=box,
                   box={'center': [1, 2, 3], 'size': [24, 24, 24]})
        condition['jobs'] = [{'id': 'box01', 'status': 'prepared', 'candidate_class_rank': 1, 'candidate_original_rank': 1,
                              'prepared_directory': 'box01', 'center_A': [1, 2, 3], 'preparation_sha256': row['preparation']['sha256']}]
    else: row['failure_reason'] = 'synthetic source preparation failure'
    row['source_condition'] = jput(root, directory + '/condition.json', condition)
    if not assembly:
        row['vina_carb_outputs'] = {str(seed): f'{base}/{profile}/{case}/{arm}/vina_carb/box01/s{seed}' for seed in validator.SEEDS}
    return row


def vc_condition(root, row, seed, plan_pin, freeze_pin, experiment, *, variant='ok'):
    directory = row['vina_carb_outputs'][str(seed)]
    record = {'schema_version': 1, 'condition': row, 'seed': seed, 'configuration_sha256': plan_pin['sha256'],
              'freeze_receipt_sha256': freeze_pin['sha256'], 'driver_source_sha256': experiment['adapter_source']['sha256'],
              'atomic_accuracy_computed': False, 'experimental_ligand_coordinates_read': False,
              'status': 'prepared_input_failure', 'full_search_started': False, 'cpu_seconds': 0, 'wall_seconds': 0}
    if row['status'] != 'prepared':
        record['error'] = row['failure_reason']
        jput(root, directory + '/execution.json', record)
        return {'status': record['status'], 'record': record}, directory + '/execution.json'
    exported = put(root, directory + '/ligand.pdbqt', b'opaque aliased ligand\n')
    put(root, directory + '/receptor.pdbqt', (root / row['receptor']['path']).read_bytes())
    export_map = jput(root, directory + '/ligand.pdbqt.map.json', {'source_sha256': row['ligand']['sha256'],
                        'ligand_sha256': exported['sha256'], 'native_coordinates_accessed': False,
                        'opaque_number': validator.OpaqueFloat('0.12345678901234567890123456789')})
    values = [('receptor', 'receptor.pdbqt'), ('ligand', 'ligand.pdbqt')]
    values += [('center_' + a, v) for a, v in zip('xyz', row['box']['center'])]
    values += [('size_' + a, v) for a, v in zip('xyz', row['box']['size'])]
    values += [('cpu', 2), ('seed', seed), ('exhaustiveness', 16), ('num_modes', 20), ('energy_range', 10),
               ('chi_coeff', 1), ('chi_cutoff', 2), ('out', 'poses.pdbqt'), ('log', 'search.log')]
    put(root, directory + '/configuration.txt', ''.join(f'{k} = {v}\n' for k, v in values).encode())
    started = variant != 'unsupported_preflight'
    if started:
        for name in ['stdout.log', 'stderr.log', 'VC_log.txt', 'search.log']:
            put(root, directory + '/' + name, b'opaque native output, not a numerical parser input\n')
        put(root, directory + '/poses.pdbqt', model_blob() if variant == 'ok' else b'opaque failed partial\n')
    record.update(status=variant, full_search_started=started, error='synthetic failure' if variant != 'ok' else None,
                  executable_sha256=experiment['native_binary_sha256'], parameters=experiment['parameters'],
                  exported_ligand_sha256=exported['sha256'], receptor_sha256=row['receptor']['sha256'],
                  command=['/synthetic/vina-carb', '--config', 'configuration.txt'], exit_code=0 if variant == 'ok' else 1)
    attempt = streams(executable_sha256=experiment['native_preflight_binary_sha256'], receptor_sha256=row['receptor']['sha256'],
                      ligand_sha256=exported['sha256'], returncode=0, timed_out=False, glycosidic_links=1, fully_supported_links=1)
    if started: record['preflight'] = attempt
    else: record['native_attempt'] = attempt
    if variant == 'ok': record.update(returned_poses=2, generated_names=['C1'])
    record['artifacts'] = {p.name: pkg.digest(p.read_bytes()) for p in (root / directory).iterdir() if p.is_file()}
    jput(root, directory + '/execution.json', record)
    return {'status': variant, 'record': record, 'pose_path': directory + '/poses.pdbqt' if variant == 'ok' else None,
            'models': 2 if variant == 'ok' else 0, 'export_map': export_map}, directory + '/execution.json'


def build_main_and_assembly(root):
    config = configuration(); modern_pin = jput(root, validator.MODERN, config)
    jput(root, validator.VINA_LOCK, {'version': '1.2.7', 'assets': {'synthetic': {'sha256': 'e' * 64}}})
    rows = []
    for profile in ['strict', 'repair']:
        for case in validator.CASES:
            for arm in validator.ARMS:
                rows.append(source_row(root, profile, case, arm, len(rows) < 19, modern_pin))
    input_pin = jput(root, validator.INPUTS, {'schema_version': 1, 'conditions': rows, 'eligible_conditions': 36,
                     'source_only_index': True, 'docking_or_accuracy_outcomes_consulted': False,
                     'primary_candidate_budget': 1, 'secondary_top5_enabled': False, 'modern_vina_configuration': modern_pin})
    adapter = put(root, 'aav_benchmark/prospective_vina_carb.py', b'synthetic adapter provenance; never executed\n')
    experiment = {'schema_version': 1, 'cases': validator.CASES, 'seeds': validator.SEEDS, 'methods': validator.ARMS,
                  'receptor_profiles': ['strict', 'repair'], 'prepared_inputs': input_pin,
                  'paired_modern_vina_policy': modern_pin, 'native_binary_sha256': 'd' * 64,
                  'native_preflight_binary_sha256': 'c' * 64, 'adapter_source': adapter,
                  'parameters': {'cpu': 2, 'exhaustiveness': 16, 'num_modes': 20, 'energy_range_kcal_mol': 10, 'chi_coeff': 1, 'chi_cutoff': 2}}
    vc_pin = jput(root, validator.VC_PLAN, experiment)
    freeze_pin = jput(root, validator.VC_FREEZE, {'configuration_sha256': vc_pin['sha256'], 'frozen_before_native_execution': True})
    jput(root, validator.MAIN_FREEZE, {'configuration_sha256': modern_pin['sha256'], 'frozen_before_new_docking': True})
    progress_rows, outcomes = [], {}
    for index, row in enumerate(rows):
        modern, path, receipt = modern_condition(root, row, config, modern_pin, input_pin['sha256'], fail_one=index == 0)
        progress_rows.append({'condition_id': row['condition_id'], 'method': 'vina_1_2_7', 'status': receipt['status'], 'execution_receipt': path})
        for seed_index, seed in enumerate(validator.SEEDS):
            outcomes[(row['condition_id'], seed, validator.ALGORITHMS[0])] = modern[seed]
            variant = 'unsupported_preflight' if index == 0 and seed_index == 0 else 'failed' if index == 0 and seed_index == 1 else 'ok'
            outcome, receipt_path = vc_condition(root, row, seed, vc_pin, freeze_pin, experiment, variant=variant)
            outcomes[(row['condition_id'], seed, validator.ALGORITHMS[1])] = outcome
            progress_rows.append({'condition_id': row['condition_id'], 'seed': seed, 'method': 'vina_carb_1_0',
                                  'status': outcome['status'], 'execution_receipt': receipt_path})
    jput(root, validator.MAIN + '/paired_execution_progress.json', {'status': 'complete', 'driver_pid': 60000,
         'planned_eligible_method_seed_slots': 216, 'terminal_eligible_seed_slots': 216, 'atomic_accuracy_evaluated': False, 'records': progress_rows})
    assembly_config = jput(root, validator.ASSEMBLY_CONFIG, config); controls = []
    for profile in ['strict_monomer', 'repair_monomer']:
        for case in ['3ZWE', '5HZA']:
            for arm in validator.ARMS:
                row = source_row(root, profile, case, arm, len(controls) < 9, assembly_config, assembly=True)
                controls.append({'row': row, 'receptor_profile': profile, 'case_id': case, 'method': arm,
                                 'condition_path': row['source_condition']['path'], 'condition_sha256': row['source_condition']['sha256'],
                                 'prepared_jobs': 1 if row['status'] == 'prepared' else 0})
    assembly_pin = jput(root, validator.ASSEMBLY_INPUTS, {'schema_version': 1, 'status': 'complete', 'conditions': [{k: v for k, v in r.items() if k != 'row'} for r in controls]})
    jput(root, validator.ASSEMBLY_FREEZE, {'configuration_sha256': assembly_config['sha256']})
    jput(root, validator.ASSEMBLY_DRIVER_FREEZE, {'configuration_sha256': assembly_config['sha256'], 'prepared_conditions_sha256': assembly_pin['sha256']})
    records = []
    for control in controls:
        row = control['row']; row['source_condition'] = {k: row['source_condition'][k] for k in ['path', 'sha256']}
        _, path, receipt = modern_condition(root, row, config, assembly_config, assembly_pin['sha256'], assembly=True)
        records.append({'condition_id': row['condition_id'], 'status': receipt['status'], 'execution_receipt': path})
    jput(root, validator.ASSEMBLY + '/execution_progress.json', {'status': 'complete', 'driver_pid': 60001,
         'planned_eligible_seed_slots': 36, 'terminal_eligible_seed_slots': 36, 'atomic_or_site_accuracy_evaluated': False, 'records': records})
    return rows, outcomes


def haddock_stages(root, directory, job_id, count):
    stages = []
    for stage_name in ['0_topoaa', '1_rigidbody', '2_flexref']:
        models = []
        for number in range(1, count + 1):
            name = f'model_{number}.pdb'
            item = put(root, directory + '/run/' + stage_name + '/' + name,
                       b'opaque HADDOCK atoms and numerical payload are deliberately not parsed\n')
            models.append({'file_name': name, 'sha256': item['sha256'], 'model_id': job_id + '/' + stage_name + '/' + name,
                           'original_model_number': number, 'atom_count': 1,
                           'score': validator.OpaqueFloat('-1.234567890123456789e+9999'), 'energies': {'opaque': 'unparsed'},
                           'actual_CNS_seed': 104729 + number, 'individual_score_rank': number, 'cluster_rank': None})
        io_pin = jput(root, directory + '/run/' + stage_name + '/io.json', {'output': [{'file_name': m['file_name'], 'score': m['score']} for m in models]})
        put(root, directory + '/run/' + stage_name + '/params.cfg', b'opaque frozen parameter source\n')
        stages.append({'stage': stage_name, 'models': models, 'actual_physical_model_files': count,
                       'io_sha256': io_pin['sha256'], 'individual_score_order': [m['model_id'] for m in models]})
    return stages


def build_haddock(root, main_rows):
    helper = put(root, 'aav_benchmark/haddock_paired_cohort.py', b'synthetic frozen HADDOCK helper; never executed\n')
    wrapper = pin(root, validator.DRIVERS['haddock'])
    policy = jput(root, 'benchmarks/haddock_orchestration_recovery_policy.json', {'source_only_fixture': True})
    cns = {'path': 'tools/unavailable-CNS', 'sha256': 'b' * 64, 'bytes': 123}
    jobs, condition_rows = [], []
    for index, source in enumerate(main_rows):
        status = 'ready' if index < 13 else 'source_parameterization_gate_failed' if index < 28 else 'assembly_protocol_unsupported' if index < 34 else 'predictor_unavailable'
        condition_rows.append({'condition_id': source['condition_id'], 'status': status, 'profile': source['receptor_profile'],
                               'case_id': source['case_id'], 'method': source['method'], 'source_condition': source['source_condition']})
        base = f"{validator.HADD}/{source['receptor_profile']}/{source['case_id']}/{source['method']}"
        input_files = []
        if status == 'ready':
            for name in ['receptor.pdb', 'glycan.pdb', 'glycan_translated.sdf', 'preflight_preparation.json', 'box_air.tbl', 'box_air.json']:
                input_files.append(put(root, base + '/data/' + name, b'opaque frozen source input\n'))
        for seed in validator.SEEDS:
            job = {'job_id': source['condition_id'] + '/s' + str(seed), 'condition_id': source['condition_id'], 'seed': seed,
                   'profile': source['receptor_profile'], 'case_id': source['case_id'], 'method': source['method'],
                   'eligible_workflow_pool': True, 'status': status, 'initial_model_slots': 20, 'maximum_refined_model_slots': 20}
            if status == 'ready':
                directory = base + '/s' + str(seed)
                job.update(directory=directory, config=put(root, directory + '/run.cfg', b'opaque frozen CNS configuration\n'), input_files=input_files)
            jobs.append(job)
    plan = {'schema_version': 1, 'seeds': validator.SEEDS, 'jobs': jobs, 'conditions': condition_rows, 'driver': helper,
            'CNS_binary': cns, 'ncores': 1, 'condition_denominator': 36, 'seed_cell_denominator': 108, 'no_accuracy_evaluation': True}
    plan_pin = jput(root, validator.HADD_PLAN, plan)
    freeze_pin = jput(root, validator.HADD_FREEZE, {'status': 'frozen_before_first_HADDOCK_cohort_pose', 'execution_authorized': True,
                         'prepared_plan': plan_pin, 'frozen_at_UTC': '2026-01-01T00:00:00+00:00'})
    orphan = next(j for j in jobs if j['job_id'] == 'strict/2RDK/fpocket/box01/s104729')
    original_failure = put(root, validator.HADD + '/execution_driver.log', b'ValueError: synthetic path.relative_to(root) failure\n')
    stages = haddock_stages(root, orphan['directory'], orphan['job_id'], 20)
    proof = {'schema_version': 1, 'job_id': orphan['job_id'], 'frozen_plan': plan_pin, 'frozen_helper': helper,
             'original_execution_freeze': freeze_pin, 'frozen_config': orphan['config'], 'original_failed_driver_log': original_failure,
             'native_completion': {'elapsed_seconds': 17, 'return_code_reason': 'original return code unknown', 'elapsed_source': 'synthetic historical log'},
             'output_proof': {'stages': stages}, 'raw_artifacts': [original_failure, orphan['config']], 'process_return_code': None, 'process_CPU_seconds': None}
    proof_pin = jput(root, validator.HADD + '/orchestration_recovery/orphan_completion_proof.json', proof)
    recovery = {'schema_version': 1, 'original_plan': plan_pin, 'original_execution_freeze': freeze_pin, 'frozen_helper': helper,
                'recovery_wrapper': wrapper, 'recovery_policy': policy, 'orphan_completion_proof': proof_pin,
                'orphan_job_id': orphan['job_id'], 'existing_terminal_receipts': [], 'seed_cell_denominator': 108,
                'original_ready_jobs': 39, 'recovered_completed_jobs': 1, 'remaining_unique_native_jobs': 38, 'native_cpus_per_job': 1,
                'no_native_rerun': True, 'no_accuracy_evaluation': True}
    recovery_pin = jput(root, validator.RECOVERY_PLAN, recovery)
    recovery_freeze_pin = jput(root, validator.RECOVERY_FREEZE, {'status': 'frozen_before_HADDOCK_orchestration_resumption',
                             'execution_authorized': True, 'recovery_plan': recovery_pin, 'maximum_parallel_jobs': 6, 'native_cpus_per_job': 1})
    adoption = {'schema_version': 1, 'status': 'verified_completed_outputs_adopted_without_rerun', 'job_id': orphan['job_id'],
                'completion_proof': proof_pin, 'recovery_plan': recovery_pin, 'recovery_freeze': recovery_freeze_pin,
                'original_plan': plan_pin, 'original_driver_failure': original_failure, 'native_rerun': False,
                'original_native_return_code': None, 'original_process_CPU_seconds': None, 'native_self_reported_elapsed_seconds': 17,
                'no_accuracy_evaluation': True}
    adoption_pin = jput(root, validator.HADD + '/orchestration_recovery/adoption_receipt.json', adoption)
    jput(root, validator.HADD + '/orchestration_recovery/recovery_preparation_receipt.json', {'recovery_plan': recovery_pin, 'no_native_execution': True})
    progress_rows = []; ready_index = 0
    for job in jobs:
        if job['status'] != 'ready':
            progress_rows.append({'job_id': job['job_id'], 'status': job['status'], 'terminal': True,
                                  'eligible_workflow_pool': True, 'physical_models': 0, 'no_native_execution': True})
            continue
        directory = job['directory']; status = 'complete_outputs_recovered' if ready_index == 0 else 'native_execution_failed' if ready_index == 1 else 'orchestration_failed' if ready_index == 2 else 'completed'
        record = {'schema_version': 1, 'job_id': job['job_id'], 'terminal': True, 'eligible_workflow_pool': True,
                  'status': status, 'prepared_plan': plan_pin, 'config': job['config'], 'no_accuracy_evaluation': True}
        if status == 'orchestration_failed':
            record.update(process_return_code=None, native_retry=False, failure={'reason': 'synthetic wrapper failure'},
                          retained_partial_stages=haddock_stages(root, directory, job['job_id'], 1))
            filename = 'wrapper_failure.json'
        else:
            actual_stages = stages if ready_index == 0 else haddock_stages(root, directory, job['job_id'], 1)
            record.update(root_freeze=freeze_pin, helper_source=helper, CNS_binary=cns, ncores=1,
                          requested_initial_model_slots=20, failed_model_slots_not_replaced=True, clustering_performed=False,
                          process_return_code=1 if status == 'native_execution_failed' else 0, elapsed_seconds=17,
                          actual_initial_models=20 if ready_index == 0 else 1, actual_refined_models=20 if ready_index == 0 else 1,
                          unfilled_initial_slots=0 if ready_index == 0 else 19, final_score_metadata_complete=status != 'native_execution_failed',
                          command=['/synthetic/haddock3', 'run.cfg'], stages=actual_stages)
            if ready_index == 0:
                record.update(native_rerun=False, process_return_code=None, process_CPU_seconds=None, external_measured_wall_seconds=None,
                              recovery_receipt=adoption_pin, execution_recovery_freeze=recovery_freeze_pin, recovery_wrapper_source=wrapper)
            filename = 'execution.json'
        receipt_pin = jput(root, directory + '/' + filename, record)
        row = {**record, 'execution_receipt': receipt_pin}
        if ready_index:
            row['resource_usage'] = jput(root, directory + '/measured_resources.json', {'schema_version': 1, 'job_id': job['job_id'],
                'native_cpus_per_job': 1, 'worker_pid': 70000 + ready_index, 'original_plan': plan_pin,
                'frozen_helper': helper, 'recovery_wrapper': wrapper, 'measured_usage': {'wall_seconds': 1, 'CPU_seconds': 1}})
        progress_rows.append(row); ready_index += 1
    jput(root, validator.HADD + '/execution_progress.json', {'schema_version': 2, 'prepared_plan': plan_pin, 'root_freeze': freeze_pin,
         'recovery_plan': recovery_pin, 'execution_recovery_freeze': recovery_freeze_pin, 'seed_cell_denominator': 108,
         'terminal_cells': 108, 'all_terminal': True, 'no_accuracy_evaluation': True, 'jobs': progress_rows})
    return plan


def build_union(root, rows, outcomes):
    pool_schedule = [{'pool_id': row['condition_id'] + '/s' + str(seed)} for row in rows for seed in validator.SEEDS]
    schedule_pin = jput(root, validator.SCHEDULE, {'schema_version': 1, 'pools': pool_schedule})
    lock = jput(root, 'benchmarks/additional_method.lock.json', {'schema_version': 1, 'synthetic_method_metadata': True})
    plan = {'schema_version': 1, 'source_inputs': pin(root, validator.INPUTS), 'source_schedule': schedule_pin,
            'source_experiment': pin(root, validator.VC_PLAN), 'source_freeze': pin(root, validator.VC_FREEZE),
            'source_modern_experiment': pin(root, validator.MODERN), 'source_modern_freeze': pin(root, validator.MAIN_FREEZE),
            'method_lock': lock, 'adapter_source': pin(root, validator.DRIVERS['union']),
            'eligible_seed_pools': 108, 'eligible_source_method_slots': 216, 'deduplicate_or_filter_poses': False}
    plan_pin = jput(root, validator.UNION_PLAN, plan)
    jput(root, validator.UNION_FREEZE, {'configuration_sha256': plan_pin['sha256'], 'frozen_before_fixed_pose_scoring': True})
    pools, poses = [], []
    for row in rows:
        for seed in validator.SEEDS:
            identity = row['condition_id']; pool_id = identity + '/s' + str(seed)
            pool = {'pool_id': pool_id, 'condition_id': identity, 'case_id': row['case_id'], 'receptor_profile': row['receptor_profile'],
                    'arm': row['method'], 'seed': seed, 'eligible_workflow_pool': True,
                    'source_statuses': {}, 'source_record_sha256': {}, 'original_source_orders': {a: [] for a in validator.ALGORITHMS}}
            for algorithm_index, algorithm in enumerate(validator.ALGORITHMS):
                outcome = outcomes[(identity, seed, algorithm)]; source_sha = validator.canonical_sha256(outcome['record'])
                pool['source_statuses'][algorithm] = outcome['status']; pool['source_record_sha256'][algorithm] = source_sha
                if outcome['status'] != 'ok': continue
                raw = (root / outcome['pose_path']).read_bytes()
                for rank, model_sha in validator.model_hashes(raw).items():
                    pose_id = pool_id + f'/{algorithm}/p{rank:02d}'
                    source = {'pose_id': pose_id, 'pool_id': pool_id, 'condition_id': identity, 'complex_id': row['case_id'],
                              'receptor_profile': row['receptor_profile'], 'arm': row['method'], 'seed': seed,
                              'source_algorithm': algorithm, 'source_algorithm_rank': algorithm_index, 'model_rank': rank,
                              'candidate_class_rank': 1, 'candidate_original_rank': 1,
                              'original_method_score': validator.OpaqueFloat('-1.234567890123456789e+9999'),
                              'source_pose_file': outcome['pose_path'], 'source_pose_file_sha256': pkg.digest(raw),
                              'source_model_sha256': model_sha, 'source_record_sha256': source_sha,
                              'source_export_map': None if algorithm_index == 0 else outcome['export_map'],
                              'source_ligand_sha256': row['ligand']['sha256'], 'source_receptor_sha256': row['receptor']['sha256'],
                              'source_preparation_sha256': row['preparation']['sha256'],
                              'source_canonical_atom_map_sha256': row['canonical_atom_map']['sha256']}
                    source['pose_identity_sha256'] = validator.canonical_sha256(source)
                    poses.append(source); pool['original_source_orders'][algorithm].append(pose_id)
            pools.append(pool)
    index = {'schema_version': 1, 'configuration_sha256': plan_pin['sha256'], 'method_lock_sha256': lock['sha256'],
             'source_inputs_sha256': plan['source_inputs']['sha256'], 'eligible_seed_pools': 108, 'eligible_source_method_slots': 216,
             'terminal_source_method_slots': 216, 'all216_source_slots_terminal_and_verified': True, 'pending_source_slots': [],
             'retained_pose_count': len(poses), 'poses': poses, 'pools': pools, 'deduplicated_poses': 0,
             'native_coordinates_used': False, 'atomic_accuracy_computed': False}
    index_pin = jput(root, validator.UNION + '/source_pool.json', index)
    records = {}; success = 0
    for number, source in enumerate(poses):
        record = {'schema_version': 1, 'source': source, 'method_lock_sha256': lock['sha256'],
                  'coordinate_optimization': False, 'native_coordinates_used': False, 'native_wall_seconds': 1, 'native_cpu_seconds': 1,
                  'exported_ligand_sha256': 'a' * 64, 'export_map_sha256': 'b' * 64, 'expected_glycosidic_links': 1}
        attempt = streams(executable_sha256='c' * 64, receptor_sha256=source['source_receptor_sha256'],
                          ligand_sha256='a' * 64, returncode=0, timed_out=False)
        if number == 0: record.update(status='pose_scoring_failed', error='synthetic retained scoring failure', native_attempt=attempt)
        else:
            record.update(status='ok', score={**attempt, 'uninterpreted_energy': validator.OpaqueFloat('1.234567890123456789e+9999')})
            success += 1
        value = jput(root, validator.UNION + '/scores/' + source['pose_id'] + '.json', record)
        records[source['pose_id']] = value['sha256']
    manifest_pin = jput(root, validator.UNION + '/scores_manifest.json', {'schema_version': 1, 'configuration_sha256': plan_pin['sha256'], 'records': records})
    jput(root, validator.UNION + '/scoring_progress.json', {'schema_version': 1, 'configuration_sha256': plan_pin['sha256'],
         'source_pool_sha256': index_pin['sha256'], 'score_manifest_sha256': manifest_pin['sha256'], 'eligible_seed_pools': 108,
         'eligible_rescore_method_slots': 216, 'indexed_pose_count': len(poses), 'recorded_pose_count': len(records),
         'successful_pose_scores': success, 'failed_pose_scores': len(records) - success,
         'all216_source_slots_terminal_and_verified': True, 'ready_for_separate_accuracy_evaluation': True, 'atomic_accuracy_computed': False})
    return index


def cohort(root):
    package_plan = helpers.fixture(root); package_plan['limits']['max_raw_bytes'] = 2_000_000
    for role, path in validator.DRIVERS.items(): put(root, path, f'synthetic {role} source metadata; never executed\n'.encode())
    rows, outcomes = build_main_and_assembly(root)
    haddock = build_haddock(root, rows); union = build_union(root, rows, outcomes)
    for names in validator.WAIT_SOURCES.values():
        for name in names:
            if not (root / name).exists():
                put(root, name, b'synthetic source wait/join contract; never executed\n')
    jput(root, 'packaging_plan.json', package_plan)
    roles = []
    for index, (role, path) in enumerate(validator.DRIVERS.items()):
        pids = [60000 + index] + (list(range(70001, 70039)) if role == 'haddock' else [])
        source = pin(root, path)
        launch = jput(root, f'validation_context/{role}-launch.json', {'schema_version': 1, 'role': role,
                      'scope': 'synthetic-only launch provenance', 'actual_known_pids': pids})
        completion = jput(root, f'validation_context/{role}-completion.json', {
            'schema_version': 1, 'kind': 'operator_observed_synchronous_tool_completion', 'role': role,
            'tool_session_id': f'synthetic-tool-session-{role}', 'exit_code': 0,
            'observed_complete': True, 'synchronous_wait_observed': True,
            'exit_code_scope': 'tool_session_only_not_original_native_process',
            'original_native_process_return_code_inferred': False,
            'launch_provenance': launch, 'driver_source': source})
        wait = jput(root, f'validation_context/{role}-source-wait.json', {
            'schema_version': 1, 'role': role, 'source_pins': [pin(root, name) for name in validator.WAIT_SOURCES[role]],
            'wait_contract': validator.WAIT_CONTRACTS[role], 'all_workers_and_child_writers_joined': True,
            'no_unjoined_writer_lifetimes': True})
        roles.append({'role': role, 'driver_source': source, 'actual_known_pids': pids,
                      'historical_pid_inventory_complete': True, 'missing_historical_pids': [],
                      'pid_inventory_reason': 'complete_recorded_inventory', 'launch_provenance': launch,
                      'operator_completion': completion, 'source_wait_proof': wait})
    context = {'schema_version': 2, 'scope': 'actual_current_host_producers_only', 'repository_root': str(root),
               'hostname': socket.gethostname(), 'source_review_and_observed_completion_account_for_writer_lifetimes': True,
               'producers': roles}
    context_pin = jput(root, 'validation_context/producer_context.json', context)
    auth = {'schema_version': 1, 'status': 'approved_metadata_only_terminal_validation', 'validation_authorized': True,
            'packaging_plan': pin(root, 'packaging_plan.json'), 'validator_sha256': pkg.digest(SOURCE.read_bytes()),
            'packager_sha256': validator.PACKAGER_SHA256, 'producer_context': context_pin,
            'source_metadata_pins': [pin(root, name) for name in validator.STATIC_NAMES]}
    jput(root, 'validation_authorization.json', auth)
    return {'plan': package_plan, 'rows': rows, 'outcomes': outcomes, 'haddock': haddock, 'union': union, 'context': context}


def run(root, *, absent=True, candidates=None):
    with patch.object(validator, 'process_absent', return_value=absent), \
         patch.object(validator, 'current_benchmark_processes', return_value=[] if candidates is None else candidates), \
         patch.object(subprocess, 'run', side_effect=AssertionError('native/reader subprocess forbidden')), \
         patch.object(subprocess, 'Popen', side_effect=AssertionError('native/reader subprocess forbidden')):
        return validator.validate(root, 'packaging_plan.json', 'validation_authorization.json', allow_synthetic_fixture=True)


def test_full_synthetic_cohort_completes_without_writes_or_float_interpretation(tmp_path):
    data = cohort(tmp_path)
    before = {p.relative_to(tmp_path).as_posix(): p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}
    result = run(tmp_path)
    assert result['counts'] == {'main': 216, 'assembly': 36, 'haddock': 108, 'union_sources': 216,
                                'pools': 108, 'retained': len(data['union']['poses']), 'recorded': len(data['union']['poses'])}
    assert result['writers_stopped'] is True and result['pending'] == []
    assert result['float_number_tokens_retained_without_numeric_conversion'] is True
    assert result['accuracy_computed'] is result['numeric_pose_or_energy_payloads_interpreted'] is False
    assert result['synthetic_test_fixture'] is True
    checks = result['producer_quiescence_checks']
    assert checks['all_actual_known_pids_absent_before_and_after'] is True
    assert checks['filtered_current_benchmark_command_candidates_absent_before_and_after'] is True
    assert checks['synchronous_completion_and_source_wait_proofs_verified'] is True
    assert checks['original_native_process_return_codes_inferred_from_tool_exit'] is False
    assert {p.relative_to(tmp_path).as_posix(): p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()} == before
    nonready = [j for j in data['haddock']['jobs'] if j['status'] != 'ready']
    assert len(nonready) == 69 and all('directory' not in j and 'config' not in j for j in nonready)
    assert len([j for j in data['haddock']['jobs'] if j['status'] == 'ready']) == 39


@pytest.fixture(scope='module')
def base_cohort(tmp_path_factory):
    root = tmp_path_factory.mktemp('readiness-base').resolve()
    previous = sys.dont_write_bytecode; sys.dont_write_bytecode = True
    try: data = cohort(root)
    finally: sys.dont_write_bytecode = previous
    return root, data


@pytest.fixture
def synthetic(tmp_path, base_cohort):
    base, data = base_cohort
    shutil.copytree(base, tmp_path, dirs_exist_ok=True)
    context = validator.load_json((tmp_path / 'validation_context/producer_context.json').read_bytes())
    context['repository_root'] = str(tmp_path)
    context_pin = jput(tmp_path, 'validation_context/producer_context.json', context)
    auth = validator.load_json((tmp_path / 'validation_authorization.json').read_bytes())
    auth['producer_context'] = context_pin; auth['validator_sha256'] = pkg.digest(SOURCE.read_bytes())
    jput(tmp_path, 'validation_authorization.json', auth)
    return tmp_path, deepcopy(data)


def mutate_json(root, name, change):
    value = validator.load_json((root / name).read_bytes()); change(value); jput(root, name, value)
    return value


def resync_modern_receipt(root, row):
    runs = validator.load_json((root / row['modern_runs']).read_bytes())
    receipt_path = str(Path(row['modern_runs']).parent.parent / 'execution.json')
    return mutate_json(root, receipt_path, lambda r: r.update(runs_sha256=pin(root, row['modern_runs'])['sha256'],
                       seed_statuses={str(v['seed']): v['status'] for v in runs}))


@pytest.mark.parametrize('change', ['missing_seed', 'duplicate_seed', 'pending_seed', 'source_failure_slots', 'vc_pending', 'vc_identity', 'progress_missing'])
def test_main216_rejects_missing_duplicate_pending_and_rebound_metadata(synthetic, change):
    root, data = synthetic; row = data['rows'][0]
    if change in {'missing_seed', 'duplicate_seed', 'pending_seed'}:
        def corrupt(rows):
            if change == 'missing_seed': rows.pop()
            elif change == 'duplicate_seed': rows[-1]['seed'] = rows[0]['seed']
            else: rows[0]['status'] = 'running'
        mutate_json(root, row['modern_runs'], corrupt); resync_modern_receipt(root, row)
    elif change == 'source_failure_slots':
        source = next(r for r in data['rows'] if r['status'] != 'prepared')
        path = str(Path(source['modern_runs']).parent.parent / 'execution.json')
        mutate_json(root, path, lambda r: r.update(eligible_seed_slots=validator.SEEDS[:-1]))
    elif change.startswith('vc_'):
        path = row['vina_carb_outputs'][str(validator.SEEDS[-1])] + '/execution.json'
        if change == 'vc_pending': mutate_json(root, path, lambda r: r.update(status='running'))
        else: mutate_json(root, path, lambda r: r.update(seed=999999))
    else: mutate_json(root, validator.MAIN + '/paired_execution_progress.json', lambda p: p['records'].pop())
    with pytest.raises(ValueError): run(root)


@pytest.mark.parametrize('change', ['missing_seed', 'duplicate_condition', 'report72', 'index_binding'])
def test_assembly36_keeps_actual_control_denominator_and_index_binding(synthetic, change):
    root, _ = synthetic
    if change == 'missing_seed':
        row = {'modern_runs': validator.ASSEMBLY + '/strict_monomer/3ZWE/fpocket/vina/box01/runs.json'}
        mutate_json(root, row['modern_runs'], lambda values: values.pop()); resync_modern_receipt(root, row)
    elif change == 'index_binding':
        mutate_json(root, validator.ASSEMBLY + '/strict_monomer/3ZWE/fpocket/vina/execution.json', lambda r: r.update(prepared_conditions_sha256='0' * 64))
    else:
        def corrupt(progress):
            if change == 'duplicate_condition': progress['records'][-1] = deepcopy(progress['records'][0])
            else: progress['planned_eligible_seed_slots'] = progress['terminal_eligible_seed_slots'] = 72
        mutate_json(root, validator.ASSEMBLY + '/execution_progress.json', corrupt)
    with pytest.raises(ValueError): run(root)


@pytest.mark.parametrize('change', ['missing', 'duplicate', 'nonready_physical_receipt', 'unknown_status', 'progress_mismatch', 'recovered_return', 'recovered_CPU', 'recovered_wall'])
def test_haddock108_progress_only_failures_and_recovery_unknowns_are_preserved(synthetic, change):
    root, _ = synthetic; path = validator.HADD + '/execution_progress.json'
    progress = validator.load_json((root / path).read_bytes())
    if change == 'missing': progress['jobs'].pop()
    elif change == 'duplicate': progress['jobs'][-1] = deepcopy(progress['jobs'][0])
    elif change == 'nonready_physical_receipt':
        next(r for r in progress['jobs'] if r.get('no_native_execution'))['execution_receipt'] = progress['jobs'][0]['execution_receipt']
    elif change == 'progress_mismatch': progress['jobs'][0]['config'] = {'path': 'other', 'sha256': '0' * 64}
    else:
        row = progress['jobs'][0]; physical_path = row['execution_receipt']['path']
        record = validator.load_json((root / physical_path).read_bytes())
        if change == 'unknown_status': record['status'] = 'running'
        elif change == 'recovered_return': record['process_return_code'] = 0
        elif change == 'recovered_CPU': record['process_CPU_seconds'] = 0
        else: record['external_measured_wall_seconds'] = 0
        progress['jobs'][0] = {**record, 'execution_receipt': jput(root, physical_path, record)}
    jput(root, path, progress)
    with pytest.raises(ValueError): run(root)


def resync_union_progress(root):
    mutate_json(root, validator.UNION + '/scoring_progress.json', lambda p: p.update(
        source_pool_sha256=pin(root, validator.UNION + '/source_pool.json')['sha256'],
        score_manifest_sha256=pin(root, validator.UNION + '/scores_manifest.json')['sha256']))


@pytest.mark.parametrize('change', ['key_substitution', 'missing_pose', 'duplicate_pose', 'duplicate_pool', 'source_record_hash',
                                  'pending_score', 'typed_source_rebind', 'orphan_score', 'stale_progress'])
def test_union_checks_exact_pose_keys_models_source_hashes_and_final_canonical_set(synthetic, change):
    root, _ = synthetic; source_path = validator.UNION + '/source_pool.json'; manifest_path = validator.UNION + '/scores_manifest.json'
    index = validator.load_json((root / source_path).read_bytes())
    if change == 'key_substitution':
        def replace(manifest):
            identity = next(iter(manifest['records'])); value = manifest['records'].pop(identity)
            manifest['records']['equal_count_wrong_identity'] = value
        mutate_json(root, manifest_path, replace)
    elif change == 'missing_pose':
        index['poses'].pop(); index['retained_pose_count'] -= 1; jput(root, source_path, index)
    elif change == 'duplicate_pose':
        index['poses'][-1] = deepcopy(index['poses'][0]); jput(root, source_path, index)
    elif change == 'duplicate_pool':
        index['pools'][-1] = deepcopy(index['pools'][0]); jput(root, source_path, index)
    elif change == 'source_record_hash':
        index['pools'][0]['source_record_sha256'][validator.ALGORITHMS[0]] = '0' * 64; jput(root, source_path, index)
    elif change in {'pending_score', 'typed_source_rebind'}:
        identity = index['poses'][0]['pose_id']; record_path = validator.UNION + '/scores/' + identity + '.json'
        def corrupt(record):
            if change == 'pending_score': record['status'] = 'preparing'
            else: record['source']['model_rank'] = True  # Python equality would otherwise equate this to rank1.
        mutate_json(root, record_path, corrupt)
        mutate_json(root, manifest_path, lambda m: m['records'].update({identity: pin(root, record_path)['sha256']}))
    elif change == 'orphan_score': jput(root, validator.UNION + '/scores/orphan.json', {'opaque': 'unmanifested evidence'})
    else: mutate_json(root, validator.UNION + '/scoring_progress.json', lambda p: p.update(source_pool_sha256='0' * 64))
    if change != 'stale_progress': resync_union_progress(root)
    with pytest.raises(ValueError): run(root)


def update_context(root, change):
    path = 'validation_context/producer_context.json'
    mutate_json(root, path, change)
    mutate_json(root, 'validation_authorization.json', lambda a: a.update(producer_context=pin(root, path)))


@pytest.mark.parametrize('change', ['wrong_host', 'wrong_root', 'scope', 'unaccounted', 'missing_role', 'duplicate_pid',
                                   'bool_pid', 'driver_rebind', 'live_pid', 'old_top_claim', 'old_pid_key', 'current_writer'])
def test_current_host_producer_context_blocks_inventory_before_any_raw_scan(synthetic, change):
    root, _ = synthetic
    if change not in {'live_pid', 'current_writer'}:
        def corrupt(context):
            if change == 'wrong_host': context['hostname'] = 'different-host.invalid'
            elif change == 'wrong_root': context['repository_root'] = '/different/root'
            elif change == 'scope': context['scope'] = 'historical-portable-pids'
            elif change == 'unaccounted': context['source_review_and_observed_completion_account_for_writer_lifetimes'] = False
            elif change == 'missing_role': context['producers'].pop()
            elif change == 'duplicate_pid': context['producers'][-1]['actual_known_pids'] = context['producers'][0]['actual_known_pids']
            elif change == 'bool_pid': context['producers'][0]['actual_known_pids'] = [True]
            elif change == 'old_top_claim': context['all_writer_descendants_accounted_for'] = True
            elif change == 'old_pid_key': context['producers'][0]['pids'] = context['producers'][0]['actual_known_pids']
            else: context['producers'][0]['driver_source']['sha256'] = '0' * 64
        update_context(root, corrupt)
    with patch.object(validator, 'packager', return_value=pkg), \
         patch.object(pkg, 'inventory', side_effect=AssertionError('raw inventory before stopped-writer gate')):
        candidates = [{'pid': 80000, 'tags': ['main_driver'], 'command_sha256': '1' * 64}] if change == 'current_writer' else []
        with pytest.raises(ValueError): run(root, absent=change != 'live_pid', candidates=candidates)


def test_missing_historical_union_pids_remain_null_with_source_wait_and_completion_proof(synthetic):
    root, _ = synthetic
    def incomplete(context):
        record = next(r for r in context['producers'] if r['role'] == 'union')
        record.update(actual_known_pids=None, historical_pid_inventory_complete=False,
                      missing_historical_pids=None, pid_inventory_reason='not_recorded_during_execution')
    update_context(root, incomplete)
    result = run(root)
    checks = result['producer_quiescence_checks']
    assert checks['actual_known_pids']['union'] is None
    assert checks['historical_pid_inventory_complete']['union'] is False
    assert checks['missing_historical_pids']['union'] is None
    assert checks['process_scan_is_not_exhaustive_historical_PID_discovery'] is True
    assert checks['original_native_process_return_codes_inferred_from_tool_exit'] is False
    assert checks['operator_completion_receipts']['union'] == pin(root, 'validation_context/union-completion.json')
    assert checks['source_wait_proofs']['union'] == pin(root, 'validation_context/union-source-wait.json')
    assert result['counts']['main'] == 216 and result['counts']['assembly'] == 36
    assert result['counts']['haddock'] == 108 and result['counts']['pools'] == 108


@pytest.mark.parametrize('change', ['null_main', 'null_assembly', 'complete_null', 'complete_empty', 'incomplete_list_missing',
                                   'incomplete_reason', 'complete_reason', 'missing_known_key', 'nonboolean_complete'])
def test_historical_inventory_unknowns_are_explicit_and_main_driver_pids_remain_required(synthetic, change):
    root, _ = synthetic
    def corrupt(context):
        role = 'main' if change == 'null_main' else 'assembly' if change == 'null_assembly' else 'union'
        record = next(r for r in context['producers'] if r['role'] == role)
        if change.startswith('null_'):
            record.update(actual_known_pids=None, historical_pid_inventory_complete=False,
                          missing_historical_pids=None, pid_inventory_reason='not_recorded_during_execution')
        elif change == 'complete_null': record['actual_known_pids'] = None
        elif change == 'complete_empty': record['actual_known_pids'] = []
        elif change == 'incomplete_list_missing': record.update(historical_pid_inventory_complete=False,
                                                               pid_inventory_reason='not_recorded_during_execution')
        elif change == 'incomplete_reason': record.update(historical_pid_inventory_complete=False, missing_historical_pids=None)
        elif change == 'complete_reason': record['pid_inventory_reason'] = 'not_recorded_during_execution'
        elif change == 'missing_known_key': record.pop('actual_known_pids')
        else: record['historical_pid_inventory_complete'] = 1
    update_context(root, corrupt)
    with patch.object(validator, 'packager', return_value=pkg), \
         patch.object(pkg, 'inventory', side_effect=AssertionError('invalid context reached raw inventory')):
        with pytest.raises(ValueError): run(root)


def update_role_proof(root, field, change, *, role='union'):
    context = validator.load_json((root / 'validation_context/producer_context.json').read_bytes())
    record = next(r for r in context['producers'] if r['role'] == role)
    name = record[field]['path']; mutate_json(root, name, change)
    update_context(root, lambda c: next(r for r in c['producers'] if r['role'] == role).update({field: pin(root, name)}))


@pytest.mark.parametrize('change', ['nonzero_exit', 'boolean_exit', 'blank_session', 'numeric_session', 'async_wait',
                                   'not_observed', 'native_exit_inferred', 'native_exit_scope', 'wrong_role',
                                   'launch_rebound', 'driver_rebound'])
def test_operator_completion_is_synchronous_tool_exit_proof_without_native_exit_inference(synthetic, change):
    root, _ = synthetic
    def corrupt(proof):
        if change == 'nonzero_exit': proof['exit_code'] = 1
        elif change == 'boolean_exit': proof['exit_code'] = False
        elif change == 'blank_session': proof['tool_session_id'] = ' '
        elif change == 'numeric_session': proof['tool_session_id'] = 123
        elif change == 'async_wait': proof['synchronous_wait_observed'] = False
        elif change == 'not_observed': proof['observed_complete'] = False
        elif change == 'native_exit_inferred': proof['original_native_process_return_code_inferred'] = True
        elif change == 'native_exit_scope': proof['exit_code_scope'] = 'original_native_process'
        elif change == 'wrong_role': proof['role'] = 'main'
        elif change == 'launch_rebound': proof['launch_provenance']['sha256'] = '0' * 64
        else: proof['driver_source']['sha256'] = '0' * 64
    update_role_proof(root, 'operator_completion', corrupt)
    with patch.object(validator, 'packager', return_value=pkg), \
         patch.object(pkg, 'inventory', side_effect=AssertionError('invalid completion reached raw inventory')):
        with pytest.raises(ValueError): run(root)


@pytest.mark.parametrize('change', ['unjoined_workers', 'unjoined_lifetimes', 'unknown_source_sha', 'missing_source',
                                   'source_order', 'wrong_wait_contract', 'wrong_role'])
def test_source_wait_proof_binds_all_exact_frozen_join_sources(synthetic, change):
    root, _ = synthetic
    def corrupt(proof):
        if change == 'unjoined_workers': proof['all_workers_and_child_writers_joined'] = False
        elif change == 'unjoined_lifetimes': proof['no_unjoined_writer_lifetimes'] = False
        elif change == 'unknown_source_sha': proof['source_pins'][0]['sha256'] = '0' * 64
        elif change == 'missing_source': proof['source_pins'].pop()
        elif change == 'source_order': proof['source_pins'].reverse()
        elif change == 'wrong_wait_contract': proof['wait_contract'] = 'unknown_wait_contract'
        else: proof['role'] = 'main'
    update_role_proof(root, 'source_wait_proof', corrupt)
    with patch.object(validator, 'packager', return_value=pkg), \
         patch.object(pkg, 'inventory', side_effect=AssertionError('invalid wait proof reached raw inventory')):
        with pytest.raises(ValueError): run(root)


@pytest.mark.parametrize('change', ['unauthorized', 'validator_pin', 'packager_pin', 'plan_pin', 'source_override', 'launch_pin'])
def test_root_authorization_and_pinned_provenance_are_required(synthetic, change):
    root, _ = synthetic
    if change == 'launch_pin':
        path = 'validation_context/main-launch.json'; (root / path).write_bytes(b'changed launch evidence')
    else:
        def corrupt(auth):
            if change == 'unauthorized': auth['validation_authorized'] = False
            elif change == 'validator_pin': auth['validator_sha256'] = '0' * 64
            elif change == 'packager_pin': auth['packager_sha256'] = '0' * 64
            elif change == 'plan_pin': auth['packaging_plan']['sha256'] = '0' * 64
            else: auth['source_metadata_pins'].pop()
        mutate_json(root, 'validation_authorization.json', corrupt)
    with pytest.raises(ValueError): run(root)


def test_opaque_canonical_hash_matches_frozen_python_number_spelling():
    ordinary = {'float': 0.10000000000000002, 'negative_zero': -0.0, 'small': 1e-06, 'nested': [1, True, None, 'text']}
    raw = json.dumps(ordinary, sort_keys=True, indent=2, allow_nan=False).encode()
    opaque = validator.load_json(raw)
    assert type(opaque['float']) is validator.OpaqueFloat
    assert validator.canonical_sha256(opaque) == pkg.digest(json.dumps(ordinary, sort_keys=True, separators=(',', ':'), allow_nan=False).encode())
    extreme = validator.load_json(b'{"uninterpreted":1.234567890123456789e+9999}')
    assert validator.canonical_bytes(extreme) == b'{"uninterpreted":1.234567890123456789e+9999}'
    with pytest.raises(ValueError): validator.canonical_bytes({'forbidden_interpreted_float': 0.5})


@pytest.mark.parametrize('raw', [b'{"a":1,"a":2}', b'{"unknown":NaN}', b'{"unknown":Infinity}', b'{"unknown":-Infinity}'])
def test_duplicate_or_nonfinite_json_is_never_tolerated(raw):
    with pytest.raises(ValueError): validator.load_json(raw)


def test_models_are_bound_as_opaque_blocks_without_atom_or_energy_parsing():
    raw = model_blob()
    assert set(validator.model_hashes(raw)) == {1, 2}
    for broken in [raw.replace(b'MODEL 2', b'MODEL 3'), raw.replace(b'ENDMDL\n', b'', 1), b'ATOM outside any model\n']:
        with pytest.raises(ValueError): validator.model_hashes(broken)


@pytest.mark.parametrize('change', ['included', 'excluded'])
def test_complete_opaque_snapshot_detects_unconsumed_source_changes(synthetic, change):
    root, _ = synthetic
    if change == 'excluded':
        plan = validator.load_json((root / 'packaging_plan.json').read_bytes())
        plan['exclusions'] = [{'path': validator.MAIN + '/global_reports', 'scope': 'tree', 'kind': 'derived', 'reason': 'Synthetic derived subtree'}]
        jput(root, 'packaging_plan.json', plan)
        mutate_json(root, 'validation_authorization.json', lambda a: a.update(packaging_plan=pin(root, 'packaging_plan.json')))
    original = validator.main_workflow
    def altered(evidence):
        value = original(evidence)
        name = validator.MAIN + ('/history/partial.bin' if change == 'included' else '/global_reports/unexpected-derived.json')
        put(root, name, b'changed opaque bytes not read by workflow projectors')
        return value
    with patch.object(validator, 'main_workflow', altered), pytest.raises(ValueError, match='snapshot changed'):
        run(root)


def test_consumed_source_metadata_is_rechecked_at_final_boundary(synthetic):
    root, _ = synthetic; original = validator.Evidence.recheck
    def changed(evidence):
        put(root, validator.MODERN, b'changed pinned source metadata')
        return original(evidence)
    with patch.object(validator.Evidence, 'recheck', changed), pytest.raises(ValueError, match='Checksum changed'):
        run(root)


def test_producer_pid_absence_is_required_again_after_validation(synthetic):
    root, _ = synthetic; calls = {}
    def absent(pid):
        calls[pid] = calls.get(pid, 0) + 1
        return calls[pid] == 1
    with patch.object(validator, 'process_absent', absent), \
         patch.object(validator, 'current_benchmark_processes', return_value=[]), \
         patch.object(subprocess, 'run', side_effect=AssertionError('native call')), \
         patch.object(subprocess, 'Popen', side_effect=AssertionError('native call')):
        with pytest.raises(ValueError, match='PID remains present'):
            validator.validate(root, 'packaging_plan.json', 'validation_authorization.json', allow_synthetic_fixture=True)
    assert any(n == 2 for n in calls.values())


def test_current_writer_command_absence_is_required_again_after_validation(synthetic):
    root, _ = synthetic
    candidate = {'pid': 81000, 'tags': ['union_driver'], 'command_sha256': '2' * 64}
    with patch.object(validator, 'process_absent', return_value=True), \
         patch.object(validator, 'current_benchmark_processes', side_effect=[[], [candidate]]) as scan, \
         patch.object(subprocess, 'run', side_effect=AssertionError('unmocked subprocess forbidden')), \
         patch.object(subprocess, 'Popen', side_effect=AssertionError('native call')):
        with pytest.raises(ValueError, match='writer command candidate'):
            validator.validate(root, 'packaging_plan.json', 'validation_authorization.json', allow_synthetic_fixture=True)
    assert scan.call_count == 2


@pytest.mark.parametrize('field', ['operator_completion', 'source_wait_proof', 'launch_provenance'])
def test_missing_role_proof_never_reaches_output_inventory(synthetic, field):
    root, _ = synthetic
    update_context(root, lambda c: c['producers'][-1].pop(field))
    with patch.object(validator, 'packager', return_value=pkg), \
         patch.object(pkg, 'inventory', side_effect=AssertionError('missing proof reached raw inventory')):
        with pytest.raises(KeyError, match=field): run(root)


@pytest.mark.parametrize('command,tag', [
    ('python3 results/structural_pose_v3/paired_execution_driver.py', 'main_driver'),
    ('python3 results/structural_assembly_pose_v3/execution_driver.py', 'assembly_driver'),
    ('python3 -m aav_benchmark.haddock_cohort_resume', 'haddock_driver'),
    ('python3 -c "score_available_union(root)"', 'union_driver'),
    ("bash -c 'python3 results/structural_pose_v3/paired_execution_driver.py'", 'main_driver'),
    ('python3 "results/structural_pose_v3/paired_execution_driver.py', 'main_driver'),
    ('cns -s unknown', 'conservative_native_writer_candidate'),
    ('fixed_pose_score --input opaque', 'conservative_native_writer_candidate'),
    ('/synthetic/root/tools/vina-1.2.7/vina --config configuration.txt', 'benchmark_native_path'),
])
def test_command_matcher_recognizes_conservative_writer_candidates_without_running_them(command, tag):
    assert tag in validator.benchmark_command_tags(Path('/synthetic/root'), command)


@pytest.mark.parametrize('command', ['', 'python3 unrelated_task.py --secret do-not-retain',
                                   'rg paired_execution_driver.py private/path',
                                   'cat results/structural_pose_v3/paired_execution_driver.py',
                                   'git diff results/haddock_paired_cohort/execution_progress.json'])
def test_command_matcher_ignores_unrelated_or_non_runner_processes(command):
    assert validator.benchmark_command_tags(Path('/synthetic/root'), command) == []


def test_ps_filter_is_mocked_sorted_and_retains_only_candidate_hashes():
    root = Path('/synthetic/root')
    secret = 'unrelated-secret-never-retained'
    command_a = 'python3 -m aav_benchmark.paired_common_pose_rescoring'
    command_b = 'cns -s opaque-native'
    self_command = 'python3 results/structural_pose_v3/paired_execution_driver.py'
    stdout = f'9002 {command_b}\n9001 {command_a}\n9003 python3 unrelated.py --secret {secret}\n{validator.os.getpid()} {self_command}\n'
    result = subprocess.CompletedProcess(['/bin/ps'], 0, stdout=stdout, stderr='opaque diagnostic')
    with patch.object(subprocess, 'run', return_value=result) as invocation, \
         patch.object(subprocess, 'Popen', side_effect=AssertionError('no actual process scan')):
        candidates = validator.current_benchmark_processes(root)
    invocation.assert_called_once_with(['/bin/ps', '-ww', '-axo', 'pid=,command='], capture_output=True, text=True, timeout=10, check=False)
    assert [row['pid'] for row in candidates] == [9001, 9002]
    assert all(set(row) == {'pid', 'tags', 'command_sha256'} for row in candidates)
    assert [row['command_sha256'] for row in candidates] == [pkg.digest(command_a.encode()), pkg.digest(command_b.encode())]
    retained = json.dumps(candidates)
    assert secret not in retained and command_a not in retained and command_b not in retained and self_command not in retained


@pytest.mark.parametrize('failure', ['nonzero', 'os_error', 'timeout', 'malformed_pid', 'missing_command'])
def test_mocked_process_inventory_errors_fail_closed(failure):
    result = subprocess.CompletedProcess(['/bin/ps'], 1 if failure == 'nonzero' else 0,
                                         stdout='not-a-pid python3 anything' if failure == 'malformed_pid' else '9001' if failure == 'missing_command' else '', stderr='')
    effect = OSError('synthetic-only ps error') if failure == 'os_error' else subprocess.TimeoutExpired(['/bin/ps'], 10) if failure == 'timeout' else None
    with patch.object(subprocess, 'run', side_effect=effect, return_value=result), \
         patch.object(subprocess, 'Popen', side_effect=AssertionError('no actual process scan')):
        with pytest.raises(ValueError): validator.current_benchmark_processes(Path('/synthetic/root'))


def test_synthetic_api_permission_is_not_available_through_cli(synthetic):
    root, _ = synthetic; output = root / 'readiness.json'
    with patch.object(validator, 'process_absent', side_effect=AssertionError('CLI must reject synthetic before PID probes')):
        with pytest.raises(ValueError, match='Real frozen structural profile'):
            validator.main(['--root', str(root), '--plan', 'packaging_plan.json', '--authorization', 'validation_authorization.json', '--output', 'readiness.json'])
    assert not output.exists()


def test_explicit_root_alias_is_refused_before_metadata_or_pids(synthetic):
    root, _ = synthetic; alias = root.parent / (root.name + '-alias'); alias.symlink_to(root, target_is_directory=True)
    try:
        with patch.object(validator, 'process_absent', side_effect=AssertionError('aliased source reached PID checks')):
            with pytest.raises(ValueError, match='Supplied source root traverses a symlink'):
                validator.validate(alias, 'packaging_plan.json', 'validation_authorization.json', allow_synthetic_fixture=True)
    finally: alias.unlink()
