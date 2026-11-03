"""Isolated, accuracy-blocked global readiness gate for four report workflows.

Validation probes never run accuracy: a private sentinel replaces every loaded
source-context alias, and metric/statistics/native-execution entrypoints fail
closed. Three 216/108/72 readers AND the fixed-union gate are mandatory.
"""
from __future__ import annotations

import argparse
import ast
from contextlib import contextmanager
import hashlib
import importlib
import inspect
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
MODULE_NAME = 'aav_benchmark.global_pose_report_gate'
SOURCE_PATH = 'aav_benchmark/global_pose_report_gate.py'
READERS = {
    'main_paired': dict(module='aav_benchmark.paired_method_report', entrypoint='generate_paired_report',
        plan='benchmarks/structural_paired_method_evaluation_v2.json',
        freeze='results/structural_pose_v3/paired_evaluation_v2_freeze.json',
        fixed_counts={'eligible_conditions': 36, 'eligible_seed_pairs': 108, 'eligible_method_outcomes': 216},
        terminal_key='terminal_method_outcomes', denominator=216),
    'haddock': dict(module='aav_benchmark.haddock_pose_report', entrypoint='generate_report',
        plan='benchmarks/haddock_pose_evaluation.json', freeze='results/haddock_paired_cohort/evaluation_freeze.json',
        fixed_counts={'eligible_cells': 108, 'native_ready_cells': 39}, terminal_key='terminal_cells', denominator=108),
    'assembly_control': dict(module='aav_benchmark.assembly_control_report', entrypoint='generate_report',
        plan='benchmarks/structural_assembly_pose_evaluation.json',
        freeze='results/structural_assembly_pose_v3/evaluation_freeze.json',
        fixed_counts={'control_conditions': 12, 'comparator_conditions': 12, 'eligible_seed_outcomes': 72},
        terminal_key='terminal_seed_outcomes', denominator=72),
}
BLOCKED_NAMES = {
    'score_pool', '_score_pool', 'protein_fits', '_protein_alignment', 'site_summary', 'accuracy_summary',
    'native_metric_controls', 'fit_receptor', 'ligand_rmsd', 'contact_pairs', 'contact_metrics', 'clash_metrics',
    'summarize_poses', 'site_metrics', 'paired_statistics', 'comparison_statistics', '_statistics', 'site_statistics',
    'conditional_statistics', 'macro_summaries', 'profile_sensitivity_statistics', 'paired_family_bootstrap', 'aggregate_runs',
    'execute_job', 'execute_all', 'run_resume_pool', 'verify_execution_freeze',
}
PENDING_KEYS = {'schema_version', 'status', 'metrics_eligible', 'counts', 'pending', 'source_sha256_inventory',
                'reference_validation', 'native_provenance', 'control_receipt_binding_field'}
PURE_MODULE = 'aav_benchmark.paired_common_pose_report'
PURE_FIXED_COUNTS = {'eligible_seed_pools': 108, 'source_method_slots': 216, 'rescore_method_slots': 216}


class GlobalGateError(ValueError):
    pass


class ForbiddenProbePath(BaseException):
    pass


class _ReadyBoundary(BaseException):
    def __init__(self, counts):
        self.counts = counts


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def pin(root, path):
    path = Path(path); path = path if path.is_absolute() else Path(root) / path
    return dict(path=str(path.relative_to(root)), sha256=sha(path))


def _checked(root, record):
    path = Path(root) / record['path']
    if sha(path) != record['sha256']:
        raise GlobalGateError('Pinned guard/report source changed: ' + record['path'])
    return path


def source_contract(path, entrypoint):
    """Expose the reviewed context boundary; execute no report or input reader."""
    tree = ast.parse(Path(path).read_text())
    entry = next((n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == entrypoint), None)
    if entry is None:
        raise GlobalGateError('Missing frozen reader entrypoint')
    calls = [n for n in ast.walk(entry) if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == '_source_context']
    if len(calls) != 1:
        raise GlobalGateError('Expected one reviewed source-context call site')
    return dict(entrypoint_line=entry.lineno, first_source_context_line=calls[0].lineno, counts_local='counts')


def make_proposal(root=ROOT, *, union_plan, union_freeze):
    root = Path(root).resolve(); readers = {}
    for name, declaration in READERS.items():
        item = dict(declaration)
        plan_path = root / item['plan']; freeze_path = root / item['freeze']
        plan, freeze = json.loads(plan_path.read_text()), json.loads(freeze_path.read_text())
        if freeze.get('evaluation_plan_sha256') != sha(plan_path):
            raise GlobalGateError('Reader plan/freeze mismatch: ' + name)
        module_path = item['module'].replace('.', '/') + '.py'
        if module_path not in plan['metric_code_sha256']:
            raise GlobalGateError('Reader itself lacks a metric-source pin')
        item.update(plan=pin(root, plan_path), freeze=pin(root, freeze_path),
                    reader_source=dict(path=module_path, sha256=plan['metric_code_sha256'][module_path]),
                    source_contract=source_contract(root / module_path, item['entrypoint']),
                    metric_sources=plan['metric_code_sha256'])
        for path, digest in item['metric_sources'].items():
            _checked(root, dict(path=path, sha256=digest))
        readers[name] = item
    union_plan_path = root / union_plan; union_freeze_path = root / union_freeze
    union = json.loads(union_plan_path.read_text()); union_receipt = json.loads(union_freeze_path.read_text())
    if (union_receipt.get('evaluation_plan_sha256') != sha(union_plan_path)
            or union.get('frozen_before_union_pose_accuracy') is not True
            or union_receipt.get('frozen_before_union_pose_accuracy') is not True):
        raise GlobalGateError('Fixed-union plan lacks its pre-accuracy freeze')
    pure_path = PURE_MODULE.replace('.', '/') + '.py'
    pure = dict(module=PURE_MODULE, entrypoint='probe_workflow', real_entrypoint='generate_report',
                plan=pin(root, union_plan_path), freeze=pin(root, union_freeze_path),
                reader_source=dict(path=pure_path, sha256=union['metric_code_sha256'][pure_path]),
                metric_sources=union['metric_code_sha256'], fixed_counts=PURE_FIXED_COUNTS)
    for path, digest in pure['metric_sources'].items():
        _checked(root, dict(path=path, sha256=digest))
    return dict(schema_version=1, status='proposal_waiting_for_independent_review_and_root_freeze',
                guard_source=pin(root, Path(__file__)), readers=readers, pure_union=pure, all_four_required=True,
                probe_isolation='one subprocess per reader; no shared monkeypatch state',
                default_mode='validation_only_no_accuracy', real_reports_require_all_four_fresh_ready=True,
                denominators_are_views_not_independent_samples=True,
                blocked_names=sorted(BLOCKED_NAMES), probe_timeout_seconds=600)


def validate_contract(root, contract_path, freeze_path):
    root = Path(root).resolve()
    contract_path = Path(contract_path); contract_path = contract_path if contract_path.is_absolute() else root / contract_path
    freeze_path = Path(freeze_path); freeze_path = freeze_path if freeze_path.is_absolute() else root / freeze_path
    authorization = json.loads(freeze_path.read_text())
    if authorization.get('status') != 'frozen_before_global_report_readiness_probes' or authorization.get('validation_only_probes_authorized') is not True:
        raise GlobalGateError('Independent review and explicit root probe freeze required')
    if authorization.get('contract', {}).get('sha256') != sha(contract_path):
        raise GlobalGateError('Global gate contract changed')
    contract = json.loads(contract_path.read_text())
    if contract['guard_source']['path'] != SOURCE_PATH or contract['guard_source']['sha256'] != sha(Path(__file__)):
        raise GlobalGateError('The running global guard is not the frozen guard source')
    _checked(root, contract['guard_source'])
    if set(contract['readers']) != set(READERS) or contract.get('all_four_required') is not True or contract.get('blocked_names') != sorted(BLOCKED_NAMES):
        raise GlobalGateError('Three frozen readers AND the pure fixed-union gate are mandatory')
    for name, spec in contract['readers'].items():
        declared = READERS[name]
        for key in ['module', 'entrypoint', 'fixed_counts', 'terminal_key', 'denominator']:
            if spec[key] != declared[key]:
                raise GlobalGateError('Reader role/count contract changed: ' + name)
        if spec['plan']['path'] != declared['plan'] or spec['freeze']['path'] != declared['freeze']:
            raise GlobalGateError('Actual frozen reader plan paths changed')
        plan = json.loads(_checked(root, spec['plan']).read_text())
        receipt = json.loads(_checked(root, spec['freeze']).read_text())
        if (receipt.get('evaluation_plan_sha256') != spec['plan']['sha256'] or
                receipt.get('frozen_before_test_pose_accuracy') is not True or receipt.get('new_test_pose_accuracy_consulted') is not False):
            raise GlobalGateError('Reader pre-accuracy freeze changed')
        if plan['metric_code_sha256'] != spec['metric_sources']:
            raise GlobalGateError('Declared transitive metric sources changed')
        for path, digest in spec['metric_sources'].items():
            _checked(root, dict(path=path, sha256=digest))
        module_path = spec['module'].replace('.', '/') + '.py'
        if spec['reader_source'] != dict(path=module_path, sha256=plan['metric_code_sha256'][module_path]):
            raise GlobalGateError('Reader source pin differs from its frozen metric plan')
        source = _checked(root, spec['reader_source'])
        if source_contract(source, spec['entrypoint']) != spec['source_contract']:
            raise GlobalGateError('Reviewed source-context boundary changed')
    pure = contract['pure_union']
    if (pure.get('module') != PURE_MODULE or pure.get('entrypoint') != 'probe_workflow' or
            pure.get('real_entrypoint') != 'generate_report' or pure.get('fixed_counts') != PURE_FIXED_COUNTS):
        raise GlobalGateError('Pure union gate contract changed')
    union = json.loads(_checked(root, pure['plan']).read_text())
    receipt = json.loads(_checked(root, pure['freeze']).read_text())
    if (receipt.get('evaluation_plan_sha256') != pure['plan']['sha256'] or
            receipt.get('frozen_before_test_pose_accuracy') is not True or receipt.get('new_test_pose_accuracy_consulted') is not False
            or receipt.get('frozen_before_union_pose_accuracy') is not True or union.get('frozen_before_union_pose_accuracy') is not True):
        raise GlobalGateError('Fixed-union pre-accuracy freeze changed')
    if union['metric_code_sha256'] != pure['metric_sources']:
        raise GlobalGateError('Pure union metric sources changed')
    for path, digest in pure['metric_sources'].items():
        _checked(root, dict(path=path, sha256=digest))
    pure_path = PURE_MODULE.replace('.', '/') + '.py'
    if pure['reader_source'] != dict(path=pure_path, sha256=union['metric_code_sha256'][pure_path]):
        raise GlobalGateError('Pure union source pin differs from its frozen metric plan')
    _checked(root, pure['reader_source'])
    return root, contract, authorization


def _counts(counts, spec, *, ready):
    if not isinstance(counts, dict) or any(not isinstance(v, int) or isinstance(v, bool) or v < 0 for v in counts.values()):
        raise GlobalGateError('Readiness counters must be nonnegative integers')
    if any(counts.get(k) != v for k, v in spec['fixed_counts'].items()):
        raise GlobalGateError('Fixed reader denominator changed')
    terminal = counts.get(spec['terminal_key'])
    if not isinstance(terminal, int) or isinstance(terminal, bool):
        raise GlobalGateError('Missing reader terminal count')
    if (ready and terminal != spec['denominator']) or (not ready and not 0 <= terminal < spec['denominator']):
        raise GlobalGateError('Reader terminal count contradicts readiness state')
    return dict(spec['fixed_counts'], **{spec['terminal_key']: terminal})


@contextmanager
def accuracy_blocked(entry, spec, extra_namespaces=(), *, allow_context_ready=True):
    """Patch aliases by function identity as well as name; restore on every exit."""
    for name in ['aav_benchmark.qualified_pose', 'aav_benchmark.metrics', 'aav_benchmark.statistics', 'aav_benchmark.heldout_pose_report']:
        importlib.import_module(name)
    namespaces = [entry.__globals__] + [vars(m) for name, m in list(sys.modules.items())
                 if name.startswith('aav_benchmark.') and m is not None] + list(extra_namespaces)
    namespaces = list({id(n): n for n in namespaces}.values())
    identities = {id(getattr(subprocess, n)) for n in ['run', 'Popen', 'call', 'check_call', 'check_output']}
    identities.add(id(os.system))
    contexts = set()
    for namespace in namespaces:
        for name, value in list(namespace.items()):
            if name == '_source_context' and callable(value):
                contexts.add(id(value))
            if (name in BLOCKED_NAMES or (inspect.isfunction(value) and value.__module__ == 'aav_benchmark.statistics')) and callable(value):
                identities.add(id(value))
    def forbidden(*args, **kwargs):
        raise ForbiddenProbePath('Accuracy/statistics/native path attempted during validation-only probe')
    def ready(*args, **kwargs):
        # Do not inspect context arguments or construct a reference graph.
        if not allow_context_ready:
            raise ForbiddenProbePath('Pure fixed-union readiness attempted reference context')
        frame = inspect.currentframe().f_back
        while frame is not None and frame.f_code is not entry.__code__:
            frame = frame.f_back
        if frame is None:
            raise GlobalGateError('Readiness boundary called outside the pinned report entrypoint')
        counts = _counts(frame.f_locals.get('counts'), spec, ready=True)
        raise _ReadyBoundary(counts)
    changed = []
    for namespace in namespaces:
        for name, value in list(namespace.items()):
            replacement = ready if name == '_source_context' or id(value) in contexts else forbidden if name in BLOCKED_NAMES or id(value) in identities else None
            if replacement is not None:
                changed.append((namespace, name, value)); namespace[name] = replacement
    native_changed = []
    for module, names in [(subprocess, ['run', 'Popen', 'call', 'check_call', 'check_output']), (os, ['system'])]:
        for name in names:
            native_changed.append((module, name, getattr(module, name))); setattr(module, name, forbidden)
    try:
        yield
    finally:
        for module, name, value in reversed(native_changed):
            setattr(module, name, value)
        for namespace, name, value in reversed(changed):
            namespace[name] = value


def probe_callable(entry, spec, root, plan, freeze, extra_namespaces=()):
    with accuracy_blocked(entry, spec, extra_namespaces):
        try:
            result = entry(root, plan, freeze, allow_pending=True)
        except _ReadyBoundary as boundary:
            return dict(status='locally_ready_accuracy_blocked', counts=boundary.counts,
                        source_context_not_entered=True, accuracy_computed=False)
    if not isinstance(result, dict) or set(result) - PENDING_KEYS or result.get('schema_version') != 1 or result.get('status') != 'pending_incomplete_workflow' or result.get('metrics_eligible') is not False:
        raise GlobalGateError('Reader returned normally without the readiness sentinel or recognized pending schema')
    pending = result.get('pending')
    if not isinstance(pending, list) or not pending or any(not isinstance(p, str) for p in pending):
        raise GlobalGateError('Pending reader must identify unfinished workflow records')
    counts = _counts(result.get('counts'), spec, ready=False)
    return dict(status='pending_incomplete_workflow', counts=counts, pending_identifier_count=len(pending),
                source_context_not_entered=True, accuracy_computed=False)


def pure_union_result(result):
    if (not isinstance(result, dict) or result.get('metrics_eligible') is not False or result.get('accuracy_computed') is not False
            or result.get('status') not in {'pending_incomplete_workflow', 'ready_for_atomic_evaluation'}
            or not isinstance(result.get('workflow_ready'), bool)):
        raise GlobalGateError('Pure union probe returned an unexpected or metric-eligible result')
    counts = result.get('counts')
    if not isinstance(counts, dict) or any(not isinstance(v, int) or isinstance(v, bool) or v < 0 for v in counts.values()):
        raise GlobalGateError('Pure union counters must be nonnegative integers')
    if any(counts.get(k) != v for k, v in PURE_FIXED_COUNTS.items()):
        raise GlobalGateError('Fixed-union denominator changed')
    ready = result['workflow_ready']
    if ready != (result['status'] == 'ready_for_atomic_evaluation'):
        raise GlobalGateError('Pure union status and ready flag disagree')
    for key in ['terminal_source_method_slots', 'terminal_rescore_method_slots']:
        if not isinstance(counts.get(key), int) or not 0 <= counts[key] <= 216:
            raise GlobalGateError('Invalid fixed-union terminal counter')
    for key in ['retained_pose_count', 'recorded_pose_count']:
        if not isinstance(counts.get(key), int) or isinstance(counts.get(key), bool) or counts[key] < 0:
            raise GlobalGateError('Fixed-union retained/recorded pose counts must be explicit')
    if ready and (counts['terminal_source_method_slots'] != 216 or counts['terminal_rescore_method_slots'] != 216
                  or counts.get('retained_pose_count') != counts.get('recorded_pose_count')):
        raise GlobalGateError('Pure union ready without all source slots and retained fixed scores')
    return dict(status='locally_ready_accuracy_blocked' if ready else 'pending_incomplete_workflow',
                counts=counts, pure_validation_only=True, source_context_not_entered=True, accuracy_computed=False)


def _import_entry(root, spec, key='entrypoint'):
    module = importlib.import_module(spec['module'])
    expected = Path(root) / (spec['module'].replace('.', '/') + '.py')
    if Path(module.__file__).resolve() != expected.resolve():
        raise GlobalGateError('Imported reader came from another repository')
    return getattr(module, spec[key])


def probe_child(root, contract_path, freeze_path, name):
    root, contract, _ = validate_contract(root, contract_path, freeze_path)
    if name not in READERS and name != 'fixed_union':
        raise GlobalGateError('Unknown reader role')
    spec = contract['pure_union'] if name == 'fixed_union' else contract['readers'][name]
    entry = _import_entry(root, spec)
    if name == 'fixed_union':
        with accuracy_blocked(entry, spec, allow_context_ready=False):
            value = entry(root, spec['plan']['path'], spec['freeze']['path'], require_freeze=True)
        result = pure_union_result(value)
    else:
        result = probe_callable(entry, spec, root, spec['plan']['path'], spec['freeze']['path'])
    return dict(reader=name, **result)


def probe_all(root, contract_path, freeze_path):
    root, contract, _ = validate_contract(root, contract_path, freeze_path)
    contract_path = Path(contract_path); contract_path = contract_path if contract_path.is_absolute() else root / contract_path
    freeze_path = Path(freeze_path); freeze_path = freeze_path if freeze_path.is_absolute() else root / freeze_path
    results = []
    for name in list(READERS) + ['fixed_union']:
        command = [sys.executable, '-m', MODULE_NAME, '--root', str(root), '--contract', str(contract_path),
                   '--freeze', str(freeze_path), '--probe-child', name]
        completed = subprocess.run(command, cwd=root, text=True, capture_output=True, check=True,
                                   timeout=contract['probe_timeout_seconds'])
        value = json.loads(completed.stdout)
        if value.get('reader') != name or value.get('accuracy_computed') is not False or value.get('source_context_not_entered') is not True or value.get('status') not in {'pending_incomplete_workflow', 'locally_ready_accuracy_blocked'}:
            raise GlobalGateError('Malformed or unblocked reader probe output')
        if name == 'fixed_union':
            if value.get('pure_validation_only') is not True:
                raise GlobalGateError('Fixed-union child was not a pure validation probe')
            pure_union_result(dict(value, status='ready_for_atomic_evaluation' if value['status'] == 'locally_ready_accuracy_blocked' else value['status'],
                                   workflow_ready=value['status'] == 'locally_ready_accuracy_blocked', metrics_eligible=False))
        else:
            _counts(value.get('counts'), contract['readers'][name], ready=value['status'] == 'locally_ready_accuracy_blocked')
        results.append(value)
    ready = all(r['status'] == 'locally_ready_accuracy_blocked' for r in results)
    return dict(schema_version=1, status='all_four_ready_accuracy_still_blocked' if ready else 'pending_global_workflow',
                all_four_ready=ready, accuracy_computed=False, real_reports_generated=False,
                contract=pin(root, contract_path), probe_freeze=pin(root, freeze_path), readers=results)


def generate_all(root, contract_path, freeze_path, output):
    root, contract, authorization = validate_contract(root, contract_path, freeze_path)
    if authorization.get('real_reports_when_all_ready_authorized') is not True:
        raise GlobalGateError('Real reports require explicit coordinator authorization in addition to readiness')
    output = Path(output); output = (output if output.is_absolute() else root / output).resolve()
    if not output.is_relative_to(root) or output.exists():
        raise GlobalGateError('Use a fresh report output directory inside the reviewed workspace')
    readiness = probe_all(root, contract_path, freeze_path)
    if not readiness['all_four_ready']:
        raise GlobalGateError('All main216 AND HADDOCK108 AND control72 AND fixed-union gates must be ready before any real accuracy')
    # Recheck immutable sources after the readiness snapshot. No probe patch exists
    # in this parent process; only this explicit mode can enter real metrics.
    validate_contract(root, contract_path, freeze_path)
    reports = {}
    for name, spec in contract['readers'].items():
        entry = _import_entry(root, spec)
        report = entry(root, spec['plan']['path'], spec['freeze']['path'], allow_pending=False)
        if report.get('status') != 'complete' or report.get('metrics_eligible') is not True:
            raise GlobalGateError('A real reader did not complete after the global gate')
        _counts(report['counts'], spec, ready=True)
        reports[name] = report
    pure = contract['pure_union']
    union = _import_entry(root, pure, 'real_entrypoint')(
        root, pure['plan']['path'], pure['freeze']['path'], main_report=reports['main_paired'], allow_pending=False)
    if union.get('status') != 'complete' or union.get('metrics_eligible') is not True:
        raise GlobalGateError('Fixed-union real report did not complete after the global gate')
    reports['fixed_union'] = union
    payloads = {name: (json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + '\n').encode()
                for name, report in reports.items()}
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix='.global-report-', dir=output.parent))
    try:
        for name, payload in payloads.items():
            with (staging / (name + '.json')).open('xb') as stream:
                stream.write(payload)
        if output.exists():
            raise GlobalGateError('Fresh report output appeared during generation')
        staging.rename(output)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return dict(readiness, status='all_four_reports_generated', accuracy_computed=True,
                real_reports_generated=True, output=str(output.relative_to(root)))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--make-proposal', type=Path)
    parser.add_argument('--union-plan', type=Path)
    parser.add_argument('--union-freeze', type=Path)
    parser.add_argument('--contract', type=Path)
    parser.add_argument('--freeze', type=Path)
    parser.add_argument('--probe-child', choices=list(READERS) + ['fixed_union'], help=argparse.SUPPRESS)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args(argv)
    if args.make_proposal:
        if not args.union_plan or not args.union_freeze:
            parser.error('Proposal requires actual frozen --union-plan and --union-freeze')
        args.make_proposal.write_text(json.dumps(make_proposal(args.root, union_plan=args.union_plan, union_freeze=args.union_freeze), indent=2, sort_keys=True, allow_nan=False) + '\n'); return 0
    if not args.contract or not args.freeze:
        parser.error('--contract and --freeze are required')
    if args.probe_child:
        result = probe_child(args.root, args.contract, args.freeze, args.probe_child)
    elif args.generate:
        if not args.output:
            parser.error('--generate requires a fresh --output directory')
        result = generate_all(args.root, args.contract, args.freeze, args.output)
    else:
        result = probe_all(args.root, args.contract, args.freeze)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)); return 0


if __name__ == '__main__':
    raise SystemExit(main())
