"""Preserving HADDOCK receipt recovery; no native cohort launch in tests."""
from copy import deepcopy
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

from aav_benchmark import haddock_cohort_resume as source_resume
from aav_benchmark.method_protocol_sources import prepare_source_fixture


class RecoveryControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix='haddock-recovery-source-')
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.root = Path(cls.temporary.name).resolve()
        prepare_source_fixture(cls.root, source_root=source_resume.frozen.ROOT, groups=('prepared', 'haddock'))
        # The wrapper records its own source path. Execute its verified private
        # copy so __file__ truthfully belongs to the private provenance root.
        source = cls.root / 'aav_benchmark/haddock_cohort_resume.py'
        spec = importlib.util.spec_from_file_location('aav_benchmark._public_resume_source_fixture', source)
        cls.resume = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.resume)
        cls.base = cls.root / 'results/haddock_paired_cohort'
        cls.plan_path = cls.base / 'prepared_plan.json'
        cls.freeze = cls.base / 'execution_freeze.json'
        cls.plan = json.loads(cls.plan_path.read_text())
        cls.job = next(j for j in cls.plan['jobs'] if j['job_id'] == cls.resume.ORPHAN_JOB)
        cls.directory = cls.root / cls.job['directory']

    def test_documented_relative_CLI_paths_become_absolute_before_frozen_call(self):
        with tempfile.TemporaryDirectory(dir=self.base) as directory:
            destination = Path(directory) / 'job'; destination.mkdir()
            relative = str(destination.relative_to(self.root))
            plan = Path(directory) / 'plan.json'
            plan.write_text(json.dumps(dict(driver=self.plan['driver'], jobs=[dict(job_id='test', directory=relative, config={})])))
            freeze = Path(directory) / 'freeze.json'; freeze.write_text('{}')
            def fake_native(plan_path, freeze_path, job_id, root):
                self.assertTrue(root.is_absolute())
                self.assertTrue(plan_path.is_absolute())
                self.assertTrue(freeze_path.is_absolute())
                result = dict(job_id=job_id, terminal=True, status='completed')
                (destination / 'execution.json').write_text(json.dumps(result))
                return result
            with patch.object(self.resume.frozen, 'execute_job', fake_native):
                result = self.resume.measured_native_worker(self.root, plan.relative_to(self.root), freeze.relative_to(self.root), 'test')
            self.assertEqual(result['result']['status'], 'completed')
            resources = json.loads((destination / 'measured_resources.json').read_text())
            self.assertGreaterEqual(resources['measured_usage']['wall_seconds'], 0)
            self.assertGreaterEqual(resources['measured_usage']['CPU_seconds'], 0)
            self.assertEqual(resources['native_cpus_per_job'], 1)

    def test_relative_pin_failure_is_reproduced_and_absolute_normalization_fixes_it(self):
        relative_plan = self.plan_path.relative_to(self.root)
        with self.assertRaises(ValueError):
            self.resume.pin(relative_plan, self.root)
        root, plan, freeze = self.resume.normalize_paths(self.root, relative_plan, self.freeze.relative_to(self.root))
        self.assertEqual(self.resume.pin(plan, root)['sha256'], self.resume.PLAN_SHA)
        self.assertEqual(self.resume.pin(freeze, root)['sha256'], self.resume.ORIGINAL_FREEZE_SHA)

    def test_existing_first_run_is_never_reexecuted_by_frozen_helper(self):
        # Exercise the actual frozen pin and existing-run guards with explicit
        # synthetic sources; public tests never invent missing installation bytes.
        with tempfile.TemporaryDirectory(prefix='haddock-existing-run-') as directory:
            root = Path(directory).resolve()
            source = root / 'synthetic-source-pin.txt'
            source.write_text('Synthetic validation fixture; no native installation.\n')
            pin = dict(path=source.name, sha256=hashlib.sha256(source.read_bytes()).hexdigest())
            job = dict(job_id='synthetic/existing', status='ready', directory='existing-run', config=pin, input_files=[pin])
            jobs = [job] + [dict(job_id=f'synthetic/unavailable/{i}', status='predictor_unavailable') for i in range(107)]
            plan = dict(ncores=1, seed_cell_denominator=108, jobs=jobs, runtime_sources=[pin],
                        **{key: pin for key in ['driver', 'proposal', 'assembly_protocol_policy', 'source_manifest',
                           'preflight_receipt', 'release_lock', 'CNS_binary', 'release_distribution', 'native_environment_freeze']})
            plan_path = root / 'synthetic-plan.json'; plan_path.write_text(json.dumps(plan))
            freeze = root / 'synthetic-freeze.json'
            freeze.write_text(json.dumps(dict(status='frozen_before_first_HADDOCK_cohort_pose', execution_authorized=True,
                                             prepared_plan=dict(sha256=hashlib.sha256(plan_path.read_bytes()).hexdigest()))))
            existing = root / job['directory']; existing.mkdir()
            for marker in ['run', 'execution.json']:
                with self.subTest(marker=marker):
                    target = existing / marker
                    if marker == 'run': target.mkdir()
                    else: target.write_text('{}')
                    with patch.object(self.resume.frozen.subprocess, 'run') as native:
                        with self.assertRaisesRegex(self.resume.frozen.PairedCohortError, 'Refuse repeated'):
                            self.resume.frozen.execute_job(plan_path, freeze, job['job_id'], root=root)
                        native.assert_not_called()
                    if target.is_dir(): target.rmdir()
                    else: target.unlink()

    @unittest.skipUnless(os.environ.get('AAV_LOCAL_NATIVE_VALIDATION') == '1',
                         'Set AAV_LOCAL_NATIVE_VALIDATION=1 for separate local release/CNS/runtime pin checks')
    def test_local_full_installation_satisfies_the_original_execution_freeze(self):
        root = source_resume.frozen.ROOT
        plan = root / 'results/haddock_paired_cohort/prepared_plan.json'
        freeze = root / 'results/haddock_paired_cohort/execution_freeze.json'
        checked = source_resume.frozen.verify_execution_freeze(plan, freeze, root=root)
        self.assertEqual(len(checked['jobs']), 108)

    def test_actual_first_outputs_are_complete_without_inventing_process_telemetry(self):
        log = self.resume.completion_log_proof(self.directory / 'stdout.log', self.freeze)
        proof = self.resume.completed_output_proof(self.directory, self.resume.ORPHAN_JOB, self.job)
        self.assertEqual(proof['actual_initial_models'], 20)
        self.assertEqual(proof['actual_refined_models'], 20)
        self.assertEqual(proof['full_model_atom_count'], 1083)
        self.assertEqual(len(proof['refinement_lineage']), 20)
        self.assertEqual(log['elapsed_seconds'], 373)
        self.assertIsNone(log['process_return_code'])
        self.assertIsNone(log['process_CPU_seconds'])
        self.assertIsNone(log['measured_external_elapsed_seconds'])

    def test_sanitized_PDB_serialization_preserves_all_physical_input_atoms(self):
        for molecule, count in [('receptor', 760), ('glycan', 23)]:
            proof = self.resume.coordinate_identity(self.directory.parent / 'data' / (molecule + '.pdb'),
                                               self.directory / 'run/data/0_topoaa' / (molecule + '.pdb'))
            self.assertEqual(proof['atom_count'], count)
            self.assertEqual(proof['maximum_coordinate_change_A'], 0)
            self.assertFalse(proof['byte_identical'])

    def _copy_outputs(self, directory):
        destination = Path(directory) / 'copy'; destination.mkdir()
        for stage in ['0_topoaa', '1_rigidbody', '2_flexref']:
            shutil.copytree(self.directory / 'run' / stage, destination / 'run' / stage)
        return destination

    def test_missing_model_and_incomplete_score_or_lineage_reject_adoption(self):
        for corruption in ['missing_model', 'score', 'lineage']:
            with tempfile.TemporaryDirectory(dir=self.base) as directory:
                destination = self._copy_outputs(directory)
                if corruption == 'missing_model':
                    next((destination / 'run/2_flexref').glob('*.pdb')).unlink()
                else:
                    path = destination / 'run/2_flexref/io.json'
                    data = json.loads(path.read_text())
                    if corruption == 'score':
                        data['output'][0]['unw_energies']['vdw'] = float('nan')
                    else:
                        data['output'][0]['ori_name'] = 'unplanned.pdb'
                    path.write_text(json.dumps(data))
                with self.assertRaises(self.resume.RecoveryError):
                    self.resume.completed_output_proof(destination, self.resume.ORPHAN_JOB, self.job)

    def test_missing_finished_marker_does_not_establish_completion(self):
        with tempfile.TemporaryDirectory(dir=self.base) as directory:
            path = Path(directory) / 'stdout'
            path.write_text((self.directory / 'stdout.log').read_text().replace('cli INFO] Finished at ', 'cli INFO] Interrupted at '))
            with self.assertRaisesRegex(self.resume.RecoveryError, 'completion'):
                self.resume.completion_log_proof(path, self.freeze)

    def test_resource_delta_excludes_other_processes_and_does_not_subtract_RSS(self):
        before = dict(worker=dict(user_seconds=2, system_seconds=1, maximum_rss=100),
                      waited_children=dict(user_seconds=10, system_seconds=5, maximum_rss=900))
        after = dict(worker=dict(user_seconds=2.5, system_seconds=1.5, maximum_rss=110),
                     waited_children=dict(user_seconds=13, system_seconds=7, maximum_rss=1000))
        measured = self.resume.measured_usage(before, after, 9)
        self.assertEqual(measured['CPU_seconds'], 6)
        self.assertEqual(measured['wall_seconds'], 9)
        self.assertEqual(measured['after']['waited_children']['maximum_rss'], 1000)
        invalid = deepcopy(after); invalid['worker']['user_seconds'] = 1
        with self.assertRaises(self.resume.RecoveryError):
            self.resume.measured_usage(before, invalid, 9)

    def test_concurrency_changes_only_after_all156_old_jobs_terminal(self):
        done = dict(scheduled_seed_runs=156, canonical_terminal_seed_runs=156, running_seed_runs=0, status='complete')
        self.assertEqual(self.resume.workflow_concurrency_limit(done, 6), 6)
        for field, value in [('scheduled_seed_runs', 155), ('canonical_terminal_seed_runs', 155),
                             ('running_seed_runs', 1), ('status', 'running')]:
            changed = dict(done); changed[field] = value
            self.assertEqual(self.resume.workflow_concurrency_limit(changed, 6), 2)
        self.assertEqual(self.resume.workflow_concurrency_limit({}, 6), 2)
        done['atomic_accuracy_evaluated'] = False
        self.assertEqual(self.resume.workflow_concurrency_limit(done, 4), 4)

    def test_control_cannot_change_native_CPUs_or_exceed_root_cap(self):
        with tempfile.TemporaryDirectory(dir=self.base) as directory:
            path = Path(directory) / 'control.json'
            path.write_text(json.dumps(dict(active_job_limit=6, stop_after_running=False)))
            self.assertEqual(self.resume.read_control(path, 6)['active_job_limit'], 6)
            with self.assertRaises(self.resume.RecoveryError):
                self.resume.read_control(path, 2)
            path.write_text(json.dumps(dict(active_job_limit=2, stop_after_running=False, native_cpus=2)))
            with self.assertRaises(self.resume.RecoveryError):
                self.resume.read_control(path, 6)

    def test_existing_terminal_receipt_is_pinned_and_changed_config_rejected(self):
        with tempfile.TemporaryDirectory(dir=self.base) as directory:
            path = Path(directory); (path / 'execution.json').write_text(json.dumps(dict(
                job_id='test', terminal=True, prepared_plan=dict(sha256=self.resume.PLAN_SHA), config=dict(sha256='expected'))))
            job = dict(job_id='test', directory=str(path.relative_to(self.root)), config=dict(sha256='expected'))
            pins = self.resume._existing_receipts(self.root, dict(jobs=[job]))
            self.assertEqual(len(pins), 1)
            contents = (path / 'execution.json').read_bytes()
            self.resume._checked(self.root / pins[0]['path'], pins[0]['sha256'])
            self.assertEqual(contents, (path / 'execution.json').read_bytes())
            job['config']['sha256'] = 'changed'
            with self.assertRaisesRegex(self.resume.RecoveryError, 'config'):
                self.resume._existing_receipts(self.root, dict(jobs=[job]))

    def test_recovery_freeze_required_before_any_adoption_or_native_call(self):
        with tempfile.TemporaryDirectory(dir=self.base) as directory:
            path = Path(directory) / 'freeze.json'; path.write_text('{}')
            with patch.object(self.resume.frozen, 'execute_job') as native:
                with self.assertRaisesRegex(self.resume.RecoveryError, 'recovery freeze'):
                    self.resume.verify_recovery_authorization(self.root, Path(directory) / 'plan.json', path)
                native.assert_not_called()

    def test_exclusive_coordinator_lock_prevents_duplicate_submission(self):
        with tempfile.TemporaryDirectory(dir=self.base) as directory:
            root = self.root.resolve(); plan = Path(directory) / 'recovery_plan.json'; plan.write_text('{}')
            lock = Path(directory) / 'coordinator.lock'; lock.write_text(json.dumps(dict(pid=999999)))
            with patch.object(self.resume, 'verify_recovery_authorization', return_value=(root, plan)), \
                 patch.object(self.resume, '_run_resume_pool_impl') as pool:
                with self.assertRaisesRegex(self.resume.RecoveryError, 'coordinator lock'):
                    self.resume.run_resume_pool(root, plan, 'unused', 'unused')
                pool.assert_not_called()
            self.assertTrue(lock.exists())


if __name__ == '__main__':
    unittest.main()
