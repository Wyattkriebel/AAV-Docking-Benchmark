"""Synthetic global gate tests; no live reader, native job or real accuracy."""
from copy import deepcopy
import json
from pathlib import Path
import subprocess
import tempfile
import types
import unittest
from unittest.mock import patch

from aav_benchmark import global_pose_report_gate as gate


def fake_entry(body, counts):
    ns = {'_source_context': lambda *a, **k: (_ for _ in ()).throw(AssertionError('real context')),
          'score_pool': lambda *a, **k: (_ for _ in ()).throw(AssertionError('real metric'))}
    code = 'def entry(root, plan, freeze, *, allow_pending=False):\n    counts = ' + repr(counts) + '\n'
    exec(code + '\n'.join('    ' + line for line in body.splitlines()), ns)
    return ns['entry'], ns


class GlobalGateControls(unittest.TestCase):
    def setUp(self):
        self.spec = gate.READERS['haddock']
        self.ready = dict(self.spec['fixed_counts'], terminal_cells=108)
        self.pending = dict(self.spec['fixed_counts'], terminal_cells=70)

    def test_sentinel_requires_full_caller_counts_and_restores_context(self):
        entry, ns = fake_entry('_source_context(object())', self.ready)
        original = ns['_source_context']
        result = gate.probe_callable(entry, self.spec, '.', 'p', 'f')
        self.assertEqual(result['status'], 'locally_ready_accuracy_blocked')
        self.assertFalse(result['accuracy_computed'])
        self.assertIs(ns['_source_context'], original)
        early, _ = fake_entry('_source_context(object())', self.pending)
        with self.assertRaises(gate.GlobalGateError):
            gate.probe_callable(early, self.spec, '.', 'p', 'f')

    def test_context_arguments_never_inspected(self):
        class Bomb:
            def __iter__(self): raise AssertionError('argument read')
            def __repr__(self): raise AssertionError('argument read')
        entry, ns = fake_entry('_source_context(argument)', self.ready); ns['argument'] = Bomb()
        self.assertEqual(gate.probe_callable(entry, self.spec, '.', 'p', 'f')['counts']['terminal_cells'], 108)

    def test_pending_is_counts_only_and_complete_without_sentinel_fails(self):
        body = "return {'schema_version':1,'status':'pending_incomplete_workflow','metrics_eligible':False,'counts':counts,'pending':['job']}"
        entry, _ = fake_entry(body, self.pending)
        result = gate.probe_callable(entry, self.spec, '.', 'p', 'f')
        self.assertEqual(result['pending_identifier_count'], 1)
        self.assertNotIn('pending', result)
        for counts in [self.ready, dict(self.pending, terminal_cells=True), dict(self.pending, eligible_cells=107)]:
            entry, _ = fake_entry(body, counts)
            with self.assertRaises(gate.GlobalGateError): gate.probe_callable(entry, self.spec, '.', 'p', 'f')
        for body in ["return {'status':'complete','metrics_eligible':True}",
                     "return {'schema_version':1,'status':'pending_incomplete_workflow','metrics_eligible':False,'counts':counts,'pending':[]}"]:
            entry, _ = fake_entry(body, self.pending)
            with self.assertRaises(gate.GlobalGateError): gate.probe_callable(entry, self.spec, '.', 'p', 'f')

    def test_metric_and_stats_aliases_blocked_and_restored(self):
        from aav_benchmark.statistics import aggregate_runs
        for name in ['metric_alias', 'stat_alias']:
            entry, ns = fake_entry(name + '([])', self.ready)
            ns['metric_alias'] = ns['score_pool']; ns['stat_alias'] = aggregate_runs
            original = dict(ns)
            with self.assertRaises(gate.ForbiddenProbePath): gate.probe_callable(entry, self.spec, '.', 'p', 'f')
            self.assertIs(ns['score_pool'], original['score_pool'])
            self.assertIs(ns['stat_alias'], original['stat_alias'])

    def test_native_alias_cannot_call_process(self):
        calls = []
        def fake_native(*a, **k): calls.append(a)
        entry, ns = fake_entry('native_alias([])', self.ready); ns['native_alias'] = fake_native
        with patch.object(subprocess, 'run', fake_native):
            with self.assertRaises(gate.ForbiddenProbePath): gate.probe_callable(entry, self.spec, '.', 'p', 'f')
        self.assertEqual(calls, [])

    def test_other_exception_propagates_exactly(self):
        error = ValueError('source corrupt'); entry, ns = fake_entry('raise error', self.ready); ns['error'] = error
        with self.assertRaises(ValueError) as caught: gate.probe_callable(entry, self.spec, '.', 'p', 'f')
        self.assertIs(caught.exception, error)

    def pure(self, ready=True):
        return dict(status='ready_for_atomic_evaluation' if ready else 'pending_incomplete_workflow',
                    workflow_ready=ready, metrics_eligible=False, accuracy_computed=False,
                    counts=dict(gate.PURE_FIXED_COUNTS, terminal_source_method_slots=216,
                                terminal_rescore_method_slots=216 if ready else 120,
                                retained_pose_count=40, recorded_pose_count=40 if ready else 20))

    def test_union_requires_every_source_and_recorded_pose_and_forbids_context(self):
        self.assertEqual(gate.pure_union_result(self.pure())['status'], 'locally_ready_accuracy_blocked')
        for field in ['terminal_source_method_slots', 'terminal_rescore_method_slots', 'recorded_pose_count']:
            value = self.pure(); value['counts'][field] -= 1
            with self.assertRaises(gate.GlobalGateError): gate.pure_union_result(value)
        missing = self.pure(); del missing['counts']['retained_pose_count']; del missing['counts']['recorded_pose_count']
        with self.assertRaises(gate.GlobalGateError): gate.pure_union_result(missing)
        entry, _ = fake_entry('_source_context(object())', self.ready)
        with self.assertRaises(gate.ForbiddenProbePath):
            with gate.accuracy_blocked(entry, self.spec, allow_context_ready=False): entry('.', 'p', 'f')

    def contract(self):
        readers = deepcopy(gate.READERS)
        for spec in readers.values():
            spec['plan'] = {'path':spec['plan']}; spec['freeze'] = {'path':spec['freeze']}
        return dict(readers=readers, pure_union=dict(module=gate.PURE_MODULE,
                    real_entrypoint='generate_report', plan=dict(path='u'), freeze=dict(path='uf')), probe_timeout_seconds=600)

    def pinned_contract(self, root):
        """Tiny source-only repository; no live input archive or reader is used."""
        def write(path, value):
            target = root / path; target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(value if isinstance(value, str) else json.dumps(value, sort_keys=True))
            return gate.pin(root, target)
        guard = write(gate.SOURCE_PATH, Path(gate.__file__).read_text())
        readers = {}
        for name, declared in gate.READERS.items():
            spec = deepcopy(declared); source_path = spec['module'].replace('.', '/') + '.py'
            source = write(source_path, 'def ' + spec['entrypoint'] + '():\n    counts = {}\n    _source_context()\n')
            plan = write(spec['plan'], {'metric_code_sha256':{source_path:source['sha256']}})
            receipt = write(spec['freeze'], dict(evaluation_plan_sha256=plan['sha256'], frozen_before_test_pose_accuracy=True,
                                                new_test_pose_accuracy_consulted=False))
            spec.update(plan=plan, freeze=receipt, reader_source=source, metric_sources={source_path:source['sha256']},
                        source_contract=gate.source_contract(root/source_path,spec['entrypoint']))
            readers[name] = spec
        pure_path = gate.PURE_MODULE.replace('.', '/') + '.py'; pure_source = write(pure_path, 'def probe_workflow(): pass\n')
        pure_plan = write('benchmarks/union.json',dict(frozen_before_union_pose_accuracy=True,metric_code_sha256={pure_path:pure_source['sha256']}))
        pure_freeze = write('results/union_freeze.json',dict(evaluation_plan_sha256=pure_plan['sha256'],frozen_before_union_pose_accuracy=True,
                                frozen_before_test_pose_accuracy=True,new_test_pose_accuracy_consulted=False))
        pure = dict(module=gate.PURE_MODULE,entrypoint='probe_workflow',real_entrypoint='generate_report',plan=pure_plan,
                    freeze=pure_freeze,reader_source=pure_source,metric_sources={pure_path:pure_source['sha256']},fixed_counts=gate.PURE_FIXED_COUNTS)
        contract = dict(guard_source=guard,readers=readers,pure_union=pure,all_four_required=True,blocked_names=sorted(gate.BLOCKED_NAMES))
        cpin = write('contract.json',contract)
        freeze = dict(status='frozen_before_global_report_readiness_probes',validation_only_probes_authorized=True,contract=cpin)
        write('freeze.json',freeze)
        return contract, freeze

    def test_frozen_source_and_contract_mutations_fail_closed(self):
        for mutation in ['source', 'role', 'guard', 'union_freeze']:
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as directory:
                root=Path(directory).resolve(); contract, freeze=self.pinned_contract(root)
                gate.validate_contract(root,'contract.json','freeze.json')
                if mutation == 'source':
                    path=root/contract['readers']['haddock']['reader_source']['path']; path.write_text(path.read_text()+'# changed\n')
                elif mutation == 'role': contract['readers']['haddock']['denominator']=107
                elif mutation == 'guard': contract['guard_source']=contract['readers']['haddock']['reader_source']
                else:
                    path=root/contract['pure_union']['freeze']['path']; receipt=json.loads(path.read_text())
                    receipt['frozen_before_union_pose_accuracy']=False; path.write_text(json.dumps(receipt))
                    contract['pure_union']['freeze']=gate.pin(root,path)
                (root/'contract.json').write_text(json.dumps(contract)); freeze['contract']=gate.pin(root,root/'contract.json')
                (root/'freeze.json').write_text(json.dumps(freeze))
                with self.assertRaises(gate.GlobalGateError): gate.validate_contract(root,'contract.json','freeze.json')

    def test_all_four_isolated_children_required_and_module_name_stable(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve(); contract = root / 'c'; freeze = root / 'f'; contract.write_text('{}'); freeze.write_text('{}')
            calls = []
            def child(command, **kwargs):
                name = command[-1]; calls.append(command)
                if name == 'fixed_union': value = gate.pure_union_result(self.pure(False))
                else:
                    spec = gate.READERS[name]
                    value = dict(status='locally_ready_accuracy_blocked', counts=dict(spec['fixed_counts'], **{spec['terminal_key']:spec['denominator']}), accuracy_computed=False, source_context_not_entered=True)
                return types.SimpleNamespace(stdout=json.dumps(dict(reader=name, **value)))
            with patch.object(gate, 'validate_contract', return_value=(root, self.contract(), {})), patch.object(gate.subprocess, 'run', child):
                result = gate.probe_all(root, contract, freeze)
            self.assertEqual(len(calls), 4)
            self.assertTrue(all(c[2] == gate.MODULE_NAME for c in calls))
            self.assertFalse(result['all_four_ready'])
            self.assertFalse(result['accuracy_computed'])

    def test_malformed_child_or_timeout_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve(); contract = root / 'c'; freeze = root / 'f'; contract.write_text('{}'); freeze.write_text('{}')
            bad = dict(reader='main_paired', status='locally_ready_accuracy_blocked', accuracy_computed=False,
                       source_context_not_entered=True, counts={'terminal_method_outcomes':0})
            with patch.object(gate, 'validate_contract', return_value=(root, self.contract(), {})), \
                 patch.object(gate.subprocess, 'run', return_value=types.SimpleNamespace(stdout=json.dumps(bad))):
                with self.assertRaises(gate.GlobalGateError): gate.probe_all(root, contract, freeze)
            with patch.object(gate, 'validate_contract', return_value=(root, self.contract(), {})), \
                 patch.object(gate.subprocess, 'run', side_effect=subprocess.TimeoutExpired('synthetic',1)):
                with self.assertRaises(subprocess.TimeoutExpired): gate.probe_all(root, contract, freeze)

    def test_one_pending_blocks_all_unpatched_real_readers(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            with patch.object(gate, 'validate_contract', return_value=(root, self.contract(), {'real_reports_when_all_ready_authorized':True})), \
                 patch.object(gate, 'probe_all', return_value={'all_four_ready':False}), patch.object(gate.importlib, 'import_module') as real:
                with self.assertRaisesRegex(gate.GlobalGateError, 'All main216'):
                    gate.generate_all(root,'p','f',root/'out')
                real.assert_not_called()

    def test_real_union_receives_only_fresh_in_memory_main_report(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve(); modules = {}; main = None
            for name,spec in gate.READERS.items():
                report = dict(status='complete',metrics_eligible=True,counts=dict(spec['fixed_counts'],**{spec['terminal_key']:spec['denominator']}))
                if name == 'main_paired': main = report
                modules[spec['module']] = types.SimpleNamespace(__file__=str(root/(spec['module'].replace('.','/')+'.py')),**{spec['entrypoint']:lambda *a,_r=report,**k:_r})
            seen = []
            def union(*args,main_report=None,**kwargs):
                seen.append(main_report); return dict(status='complete',metrics_eligible=True)
            modules[gate.PURE_MODULE] = types.SimpleNamespace(__file__=str(root/(gate.PURE_MODULE.replace('.','/')+'.py')),generate_report=union)
            with patch.object(gate,'validate_contract',return_value=(root,self.contract(),{'real_reports_when_all_ready_authorized':True})), \
                 patch.object(gate,'probe_all',return_value=dict(all_four_ready=True)), \
                 patch.object(gate.importlib,'import_module',side_effect=lambda name:modules[name]):
                result = gate.generate_all(root,'p','f',root/'out')
            self.assertIs(seen[0],main)
            self.assertTrue(result['real_reports_generated'])
            self.assertEqual(len(list((root/'out').glob('*.json'))),4)

    def test_real_mode_needs_explicit_authorization_and_fresh_workspace_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            with patch.object(gate,'validate_contract',return_value=(root,self.contract(),{})), patch.object(gate,'probe_all') as probe:
                with self.assertRaises(gate.GlobalGateError): gate.generate_all(root,'p','f',root/'out')
                probe.assert_not_called()

    def test_nan_rejected_before_any_partial_report_is_published(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve(); modules = {}
            for spec in gate.READERS.values():
                report = dict(status='complete',metrics_eligible=True,
                    counts=dict(spec['fixed_counts'],**{spec['terminal_key']:spec['denominator']}))
                modules[spec['module']] = types.SimpleNamespace(__file__=str(root/(spec['module'].replace('.','/')+'.py')),**{spec['entrypoint']:lambda *a,_r=report,**k:_r})
            modules[gate.PURE_MODULE] = types.SimpleNamespace(__file__=str(root/(gate.PURE_MODULE.replace('.','/')+'.py')),generate_report=lambda *a,**k:dict(status='complete',metrics_eligible=True,value=float('nan')))
            with patch.object(gate,'validate_contract',return_value=(root,self.contract(),{'real_reports_when_all_ready_authorized':True})), \
                 patch.object(gate,'probe_all',return_value=dict(all_four_ready=True)), \
                 patch.object(gate.importlib,'import_module',side_effect=lambda name:modules[name]):
                with self.assertRaises(ValueError): gate.generate_all(root,'p','f',root/'out')
            self.assertFalse((root/'out').exists())
            self.assertEqual(list(root.glob('.global-report-*')),[])
            (root/'out').mkdir(); (root/'out/already.json').write_text('{}')
            with patch.object(gate,'validate_contract',return_value=(root,self.contract(),{'real_reports_when_all_ready_authorized':True})), patch.object(gate,'probe_all') as probe:
                with self.assertRaises(gate.GlobalGateError): gate.generate_all(root,'p','f',root/'out')
                probe.assert_not_called()

    def test_actual_three_source_boundaries_exist_without_invocation(self):
        for spec in gate.READERS.values():
            path = gate.ROOT/(spec['module'].replace('.','/')+'.py')
            record = gate.source_contract(path,spec['entrypoint'])
            self.assertGreater(record['first_source_context_line'],record['entrypoint_line'])


if __name__ == '__main__':
    unittest.main()
