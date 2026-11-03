"""Actual paired receptor topology inventory and source-only amendment tests."""
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

from aav_benchmark.haddock_cohort_preparation import ROOT, _checked, pin
from aav_benchmark.method_protocol_sources import prepare_source_fixture, restore_preflight_receptor_aliases
from scripts import preflight_haddock_pairs as paired
from scripts.preflight_haddock_terminals import observed_terminal_rule


class PairedTopologyControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix='haddock-paired-topology-source-')
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.root = Path(cls.temporary.name).resolve()
        prepare_source_fixture(cls.root, source_root=ROOT, groups=('prepared', 'haddock'))
        restore_preflight_receptor_aliases(cls.root)
        cls.output = cls.root / paired.OUTPUT.relative_to(ROOT)
        cls.policy = json.loads((cls.root / paired.POLICY.relative_to(ROOT)).read_text())
        cls.amendment = json.loads((cls.output / 'observed_termini_v2/summary.json').read_text())

    def test_all_twelve_retain_glycan_and_all_source_assembly_chains(self):
        self.assertEqual(len(self.amendment['cases']), 12)
        for row in self.amendment['cases']:
            audit = row['audit']
            self.assertEqual(audit['process_return_code'], 0)
            self.assertEqual(audit['glycan']['status'], 'accepted_glycan_parameterization')
            self.assertEqual(audit['receptor']['input_chains'], audit['receptor']['output_chains'])
            self.assertTrue(audit['glycan']['post_CNS_geometry']['valid'])
            self.assertTrue(audit['no_docking_performed'])
        trimer = next(r for r in self.amendment['cases'] if r['case_id'] == '3ZWE' and r['profile'] == 'strict')
        self.assertEqual(trimer['audit']['receptor']['output_chains'], ['A', 'B', 'C'])

    def test_actual_terminal_amendment_preserves_observed_oxygens(self):
        for profile in ['strict', 'repair']:
            original = json.loads((self.output / profile / '3ZWE/paired_topology_audit.json').read_text())
            amended = next(r for r in self.amendment['cases'] if r['case_id'] == '3ZWE' and r['profile'] == profile)
            self.assertEqual(len(original['receptor']['input_only_heavy_atoms']), 3)
            self.assertTrue(all(a[-1] == 'OXT' for a in original['receptor']['input_only_heavy_atoms']))
            self.assertEqual(amended['audit']['receptor']['status'], 'exact_named_heavy_receptor_preserved')
            self.assertEqual(amended['audit']['receptor']['output_heavy_atoms'], 1983)

    def test_exact_receptor_availability_is_explicit_and_denominator_is_retained(self):
        accepted = {(r['profile'], r['case_id']) for r in self.amendment['cases']
                    if r['audit']['receptor']['status'] == 'exact_named_heavy_receptor_preserved'}
        self.assertEqual(accepted, {(p, c) for p in ['strict', 'repair'] for c in ['2RDK', '3ZWE', '4R9F']} | {('repair', '5V6F')})
        for row in self.amendment['cases']:
            if row['case_id'] in ['5HZA', '8I5U']:
                self.assertGreater(row['audit']['receptor']['shared_heavy_max_change_A'], .5)
                if row['profile'] == 'strict':
                    self.assertGreater(len(row['audit']['receptor']['CNS_added_or_renamed_heavy_atoms']), 0)

    def test_observed_terminal_rule_rejects_mixed_duplicate_or_nonterminal_inventory(self):
        atom = lambda chain, number, name: dict(chain=chain, resnum=str(number), atom_name=name)
        self.assertFalse(observed_terminal_rule([atom('A', 1, 'C'), atom('B', 1, 'C')]))
        self.assertTrue(observed_terminal_rule([atom('A', 1, 'OXT'), atom('B', 1, 'OXT')]))
        for atoms in [[atom('A', 1, 'OXT'), atom('B', 1, 'C')],
                      [atom('A', 1, 'OXT'), atom('A', 1, 'OXT')],
                      [atom('A', 1, 'OXT'), atom('A', 2, 'C')]]:
            with self.assertRaises(ValueError):
                observed_terminal_rule(atoms)

    def test_all_pre_execution_source_and_driver_pins_replay(self):
        plan = json.loads((self.output / 'execution_plan.json').read_text())
        _checked(self.root / plan['policy']['path'], plan['policy']['sha256'])
        _checked(self.root / plan['driver']['path'], plan['driver']['sha256'])
        for job in plan['jobs']:
            for key in ['receptor', 'independent_sdf', 'input_PDB', 'config', 'preparation']:
                _checked(self.root / job[key]['path'], job[key]['sha256'])
        terminal = json.loads((self.output / 'observed_termini_v2/execution_plan.json').read_text())
        _checked(self.root / terminal['policy']['path'], terminal['policy']['sha256'])
        _checked(self.root / terminal['driver']['path'], terminal['driver']['sha256'])

    def test_moved_observed_protein_atom_cannot_pass_exact_matched_input_gate(self):
        with tempfile.TemporaryDirectory(dir=self.output) as directory:
            target = Path(directory) / '2RDK'
            shutil.copytree(self.output / 'strict/2RDK', target)
            path = target / 'run-topology/0_topoaa/receptor_haddock.pdb'
            lines = path.read_text().splitlines()
            index = next(i for i, line in enumerate(lines) if line.startswith('ATOM') and line[76:78].strip() == 'C')
            line = lines[index]
            lines[index] = line[:30] + f'{float(line[30:38]) + .2:8.3f}' + line[38:]
            path.write_text('\n'.join(lines) + '\n')
            # audit_receptor has no root argument; bind only its provenance-pin
            # alias to the private fixture while preserving all numerical gates.
            with patch.object(paired, 'pin', lambda path: pin(path, self.root)):
                self.assertEqual(paired.audit_receptor(target, self.policy)['status'], 'CNS_receptor_preparation_differs')


if __name__ == '__main__':
    unittest.main()
