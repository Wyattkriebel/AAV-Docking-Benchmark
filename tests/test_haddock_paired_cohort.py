"""Prepared box/AIR cohort controls; never execute cohort docking."""
from collections import Counter
from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from aav_benchmark import haddock_paired_cohort as cohort
from aav_benchmark.method_protocol_sources import prepare_source_fixture


def synthetic_execution_plan(root, jobs):
    """Exercise frozen pin/launch guards without claiming an installation."""
    source = root / 'synthetic-source-pin.txt'
    source.write_text('Synthetic validation fixture; no executable or native installation.\n')
    record = dict(path=source.name, sha256=hashlib.sha256(source.read_bytes()).hexdigest())
    plan = dict(ncores=1, seed_cell_denominator=108, jobs=jobs,
                runtime_sources=[record], **{key: record for key in
                ['driver', 'proposal', 'assembly_protocol_policy', 'source_manifest', 'preflight_receipt',
                 'release_lock', 'CNS_binary', 'release_distribution', 'native_environment_freeze']})
    path = root / 'synthetic-prepared-plan.json'
    path.write_text(json.dumps(plan))
    freeze = root / 'synthetic-freeze.json'
    freeze.write_text(json.dumps(dict(status='frozen_before_first_HADDOCK_cohort_pose', execution_authorized=True,
                                     prepared_plan=dict(sha256=hashlib.sha256(path.read_bytes()).hexdigest()))))
    return path, freeze


class PairedCohortSourceControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix='haddock-paired-cohort-source-')
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.root = Path(cls.temporary.name).resolve()
        prepare_source_fixture(cls.root, source_root=cohort.ROOT, groups=('prepared', 'haddock'))
        cls.base = cls.root / cohort.DEFAULT_OUTPUT
        cls.plan_path = cls.base / 'prepared_plan.json'
        cls.plan = json.loads(cls.plan_path.read_text())

    def test_whole36_by3_denominator_and_method_specific_support_preserved(self):
        self.assertEqual(len(self.plan['conditions']), 36)
        self.assertEqual(len(self.plan['jobs']), 108)
        self.assertEqual(Counter(j['status'] for j in self.plan['jobs']),
                         {'ready': 39, 'source_parameterization_gate_failed': 45, 'predictor_unavailable': 6,
                          'assembly_protocol_unsupported': 18})
        self.assertTrue(all(j['eligible_workflow_pool'] for j in self.plan['jobs']))
        four_r9f = [c for c in self.plan['conditions'] if c['case_id'] == '4R9F']
        self.assertEqual(len(four_r9f), 6)
        self.assertTrue(all(c['status'] == 'ready' for c in four_r9f))
        self.assertEqual({j['seed'] for j in self.plan['jobs']}, set(cohort.SEEDS))

    def test_actual13_AIRs_parse_in_CNS_with_exact_source_atom_counts(self):
        receipt = json.loads((self.base / 'AIR_syntax_preflight/summary.json').read_text())
        self.assertTrue(receipt['no_pose_generation'])
        self.assertFalse(receipt['topology_or_coordinates_modified'])
        self.assertEqual(len(receipt['conditions']), 13)
        self.assertTrue(all(c['status'] == 'passed' and c['actual_selection_counts'] == c['expected_selection_counts']
                            for c in receipt['conditions']))
        for condition in self.plan['conditions']:
            if condition['status'] != 'ready':
                continue
            metadata = json.loads((self.root / condition['AIR_metadata']['path']).read_text())
            air = (self.root / condition['AIR']['path']).read_text()
            self.assertEqual(air.count('assign'), 1)
            self.assertEqual(air.count(' and name '), len(metadata['protein_heavy_atom_keys']) + len(metadata['full_glycan_heavy_atom_keys']))
            self.assertFalse(metadata['native_contacts_used'])
            self.assertFalse(metadata['residue_expansion'])

    def test_closed_box_union_includes_boundary_excludes_H_and_never_expands_residue(self):
        def atom(number, name, x, element='C'):
            return dict(chain='A', resnum=str(number), resname='ALA', atom_name=name,
                        element=element, xyz_A=[x, 0, 0])
        atoms = [atom(1, 'CA', 1), atom(1, 'CB', 1.001), atom(2, 'CA', -1), atom(3, 'H', 0, 'H')]
        glycan = [dict(chain='L', resnum='1', resname='MAN', atom_name='C1', element='C', xyz_A=[50, 0, 0])]
        topology = dict(atoms={str(i): dict(segment='A', resnum=a['resnum'], resname=a['resname'], atom_name=a['atom_name'], mass=12)
                               for i, a in enumerate(atoms) if a['element'] != 'H'})
        _, metadata = cohort.box_air(atoms, glycan, topology, dict(center=[0, 0, 0], size=[2, 2, 2]))
        self.assertEqual(metadata['protein_heavy_atom_keys'], [['A', '1', 'ALA', 'CA'], ['A', '2', 'ALA', 'CA']])
        with self.assertRaisesRegex(cohort.PairedCohortError, 'empty_box'):
            cohort.box_air(atoms, glycan, topology, dict(center=[10, 10, 10], size=[1, 1, 1]))
        wrong = deepcopy(topology); wrong['atoms'].clear()
        with self.assertRaisesRegex(cohort.PairedCohortError, 'CNS topology'):
            cohort.box_air(atoms, glycan, wrong, dict(center=[0, 0, 0], size=[2, 2, 2]))

    def test_original_predictor_center_survives_Vina_specific_preparation_failure(self):
        condition = dict(status='preparation_failed', jobs=[dict(status='unsupported', center_A=[99, 99, 99])],
                         candidate_selection=dict(candidates=[dict(rank=1, center_A=[1, 2, 3]), dict(rank=2, center_A=[4, 5, 6])]))
        self.assertEqual(cohort.source_box(condition, dict(box_size=[27, 27, 27]), 'fpocket')['center'], [1, 2, 3])
        self.assertEqual(cohort.source_box(condition, dict(box_size=[27, 27, 27]), 'oracle')['center'], [99, 99, 99])
        condition['candidate_selection']['candidates'] = []
        self.assertIsNone(cohort.source_box(condition, dict(box_size=27), 'p2rank'))

    def test_all39_native_configs_are_bounded_no_clustering_or_reference_modules(self):
        for job in self.plan['jobs']:
            if job['status'] != 'ready':
                continue
            config = (self.root / job['config']['path']).read_text()
            self.assertIn('ncores = 1', config)
            self.assertIn('sampling = 20', config)
            self.assertIn('sampling_factor = 1', config)
            self.assertEqual(config.count('randremoval = false'), 2)
            self.assertNotIn('caprieval', config)
            self.assertNotIn('rmsdmatrix', config)
            self.assertNotIn('clust', config)
            self.assertNotIn('reference', config)
            self.assertIn('iniseed = ' + str(job['seed']), config)
        probes = json.loads((self.base / 'setup_probes/summary.json').read_text())
        self.assertEqual(len(probes), 3)
        self.assertTrue(all(p['return_code'] == 0 and p['no_native_pose_execution'] for p in probes))

    def test_multichain_assembly_is_retained_as_adapter_failure_without_chain_removal(self):
        trimer = [c for c in self.plan['conditions'] if c['case_id'] == '3ZWE']
        self.assertEqual(len(trimer), 6)
        self.assertTrue(all(c['status'] == 'assembly_protocol_unsupported' for c in trimer))
        self.assertTrue(all(c['assembly_protocol']['source_receptor_chains'] == ['A', 'B', 'C'] for c in trimer))
        self.assertFalse(cohort.assembly_protocol_support([dict(chain='X', element='C'), dict(chain='Y', element='C')])['supported'])
        self.assertTrue(cohort.assembly_protocol_support([dict(chain='X', element='C')])['supported'])

    def test_external_exact_freeze_and_every_source_pin_required_before_execution(self):
        with tempfile.TemporaryDirectory(dir=self.base) as directory:
            freeze = Path(directory) / 'freeze.json'
            freeze.write_text(json.dumps(dict(status='not_frozen', execution_authorized=False)))
            with patch.object(cohort.subprocess, 'run') as native:
                with self.assertRaisesRegex(cohort.PairedCohortError, 'Root execution freeze'):
                    cohort.execute_job(self.plan_path, freeze, self.plan['jobs'][0]['job_id'], root=self.root)
                native.assert_not_called()
            freeze.write_text(json.dumps(dict(status='frozen_before_first_HADDOCK_cohort_pose', execution_authorized=True,
                                               prepared_plan=dict(sha256='0' * 64))))
            with self.assertRaisesRegex(ValueError, 'checksum'):
                cohort.verify_execution_freeze(self.plan_path, freeze, root=self.root)

        # The public fixture preserves all positive pin/unavailable-job controls
        # with explicit synthetic sources. No omitted installation pin is waived.
        with tempfile.TemporaryDirectory(prefix='haddock-synthetic-freeze-') as directory:
            root = Path(directory).resolve()
            jobs = [dict(job_id=f'synthetic/{i}', status='predictor_unavailable') for i in range(108)]
            plan, freeze = synthetic_execution_plan(root, jobs)
            checked = cohort.verify_execution_freeze(plan, freeze, root=root)
            self.assertEqual(len(checked['jobs']), 108)
            with patch.object(cohort.subprocess, 'run') as native:
                result = cohort.execute_job(plan, freeze, jobs[0]['job_id'], root=root)
                self.assertTrue(result['terminal'])
                self.assertEqual(result['physical_models'], 0)
                native.assert_not_called()
            (root / 'synthetic-source-pin.txt').write_text('Changed synthetic source pin.\n')
            with self.assertRaisesRegex(ValueError, 'checksum'):
                cohort.verify_execution_freeze(plan, freeze, root=root)

    @unittest.skipUnless(os.environ.get('AAV_LOCAL_NATIVE_VALIDATION') == '1',
                         'Set AAV_LOCAL_NATIVE_VALIDATION=1 for separate local release/CNS/runtime pin checks')
    def test_local_full_installation_satisfies_every_original_execution_pin(self):
        source_plan = cohort.ROOT / cohort.DEFAULT_OUTPUT / 'prepared_plan.json'
        source_freeze = cohort.ROOT / cohort.DEFAULT_OUTPUT / 'execution_freeze.json'
        checked = cohort.verify_execution_freeze(source_plan, source_freeze, root=cohort.ROOT)
        self.assertEqual(len(checked['jobs']), 108)

    def test_inventory_retains_full_precision_scores_all_models_and_stable_tie_order(self):
        with tempfile.TemporaryDirectory(dir=self.base) as directory:
            run = Path(directory); stage = run / '02_flexref'; stage.mkdir()
            outputs = []
            pdb = 'ATOM      1  C1  MAN L   1       0.000   0.000   0.000  1.00  0.00           C  \nEND\n'
            for number, score in [(10, -1.234567), (2, -1.234567), (1, -.1)]:
                name = f'flexref_{number}.pdb'; (stage / name).write_text(pdb)
                outputs.append(dict(file_name=name, score=score, seed=104729 + number, unw_energies=dict(vdw=score)))
            (stage / 'io.json').write_text(json.dumps(dict(output=outputs)))
            inventory = cohort.model_inventory(run, 'CASE/ARM/sSEED')[0]
            self.assertEqual(inventory['actual_physical_model_files'], 3)
            self.assertTrue(inventory['individual_score_order'][0].endswith('flexref_2.pdb'))
            self.assertTrue(inventory['individual_score_order'][1].endswith('flexref_10.pdb'))
            self.assertEqual(next(m for m in inventory['models'] if m['original_model_number'] == 10)['score'], -1.234567)
            self.assertTrue(all(m['cluster_rank'] is None for m in inventory['models']))


if __name__ == '__main__':
    unittest.main()
