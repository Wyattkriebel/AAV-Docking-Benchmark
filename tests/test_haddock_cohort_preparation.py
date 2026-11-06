"""Actual source-only CNS glycan topology gates and corruption controls."""
from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest

import gemmi
import numpy as np
from rdkit import Chem

from aav_benchmark.haddock_cohort_preparation import (
    ROOT, POLICY, DEFAULT_OUTPUT, CohortPreparationError, _checked,
    audit_topology, prepare_glycan,
)
from aav_benchmark.haddock_method import pdb_atoms
from aav_benchmark.method_protocol_sources import prepare_source_fixture


class ActualGlycanTopologyControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix='haddock-glycan-source-')
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.root = Path(cls.temporary.name).resolve()
        prepare_source_fixture(cls.root, source_root=ROOT, groups=('prepared', 'haddock'))
        cls.output = cls.root / DEFAULT_OUTPUT
        cls.policy = json.loads((cls.root / POLICY).read_text())

    def test_all_six_actual_CNS_topologies_preserve_complete_parents(self):
        for case in self.policy['cases']:
            with self.subTest(case=case):
                actual = audit_topology(case, self.output / case, root=self.root)
                frozen = json.loads((self.output / case / 'topology_audit.json').read_text())
                self.assertEqual(actual['status'], 'accepted_glycan_parameterization')
                self.assertEqual(actual['actual_heavy_bonds'], frozen['actual_heavy_bonds'])
                self.assertEqual(actual['input_heavy_count'], actual['typed_heavy_count'])
                self.assertEqual(actual['output_heavy_count'], actual['typed_heavy_count'])
                self.assertGreater(actual['CNS_hydrogens_built'], 0)
                self.assertEqual(actual['heavy_coordinate_max_change_A'], 0)
                self.assertLess(abs(actual['topology_total_charge']), 1e-6)
                self.assertTrue(actual['post_CNS_geometry']['valid'])
                self.assertTrue(actual['no_docking_performed'])
                self.assertEqual(frozen['process_return_code'], 0)

    def test_pre_execution_plan_pins_actual_inputs_configs_policy_and_driver(self):
        plan = json.loads((self.output / 'execution_plan.json').read_text())
        self.assertTrue(plan['written_before_first_CNS_job'])
        self.assertTrue(plan['no_docking_stages'])
        pins = [plan['policy'], plan['algorithm']]
        for case in plan['inputs']:
            pins.extend(case[key] for key in ['input_PDB', 'config', 'preparation'])
        for record in pins:
            _checked(self.root / record['path'], record['sha256'])
        # The public archive retains the executable pin as provenance metadata;
        # it does not redistribute or validate a CNS installation.
        lock = json.loads((self.root / self.policy['release_lock_path']).read_text())
        self.assertEqual(plan['CNS_binary'], lock['cns_binary'])
        for case in self.policy['cases']:
            text = (self.output / case / 'topology.cfg').read_text()
            self.assertEqual([line for line in text.splitlines() if line.startswith('[')], ['[topoaa]'])
            self.assertIn('ncores = 1', text)

    @unittest.skipUnless(os.environ.get('AAV_LOCAL_NATIVE_VALIDATION') == '1',
                         'Set AAV_LOCAL_NATIVE_VALIDATION=1 for the separate local CNS installation pin check')
    def test_local_CNS_executable_matches_the_retained_pre_execution_pin(self):
        plan = json.loads((self.output / 'execution_plan.json').read_text())
        record = plan['CNS_binary']
        _checked(ROOT / record['path'], record['sha256'])

    def test_conversion_uses_every_independent_mapped_heavy_atom_without_fitting(self):
        with tempfile.TemporaryDirectory(dir=self.output) as directory:
            prepared = prepare_glycan('3ZWE', Path(directory), root=self.root)
            molecule = Chem.SDMolSupplier(str(self.root / prepared['input_sdf']['path']), removeHs=False)[0]
            pdb = pdb_atoms(Path(directory) / 'data/glycan.pdb')
            self.assertEqual(len(pdb), molecule.GetNumHeavyAtoms())
            self.assertEqual(len(pdb), 44)
            for atom, converted in zip(prepared['atom_map'], pdb):
                self.assertEqual(atom['canonical_atom_key'], molecule.GetAtomWithIdx(atom['rdkit_index']).GetProp('canonical_atom_key'))
                self.assertEqual(converted['atom_name'], atom['atom_name'])
                self.assertEqual(converted['resname'], atom['ccd_id'])
                self.assertEqual(converted['resnum'], atom['branch_num'])
                self.assertLessEqual(np.linalg.norm(np.array(molecule.GetConformer().GetAtomPosition(atom['rdkit_index'])) - converted['xyz_A']), .000867)

    def _mutated_case(self, operation):
        directory = tempfile.TemporaryDirectory(dir=self.output)
        target = Path(directory.name) / '2RDK'
        shutil.copytree(self.output / '2RDK', target)
        stage = target / 'run-topology/0_topoaa'
        if operation in ['missing_atom', 'moved_atom']:
            path = stage / 'glycan_haddock.pdb'
            lines = path.read_text().splitlines()
            index = next(i for i, line in enumerate(lines) if line.startswith('ATOM') and line[76:78].strip() == 'C')
            if operation == 'missing_atom':
                del lines[index]
            else:
                line = lines[index]
                lines[index] = line[:30] + f'{float(line[30:38]) + .5:8.3f}' + line[38:]
            path.write_text('\n'.join(lines) + '\n')
        else:
            path = stage / 'glycan_haddock.psf'
            document = gemmi.cif.read_file(str(path))
            if operation == 'missing_bond':
                atoms = document[0].get_mmcif_category('_cns_mtf_atom.')
                heavy = {atoms['id'][i] for i, mass in enumerate(atoms['atom_mass']) if float(mass) > 2.5}
                bonds = document[0].get_mmcif_category('_cns_mtf_bond.')
                index = next(i for i, left in enumerate(bonds['id[1]']) if left in heavy and bonds['id[2]'][i] in heavy)
                document[0].set_mmcif_category('_cns_mtf_bond.', {k: [v for i, v in enumerate(values) if i != index] for k, values in bonds.items()})
            else:
                atoms = document[0].get_mmcif_category('_cns_mtf_atom.')
                index = next(i for i, mass in enumerate(atoms['atom_mass']) if float(mass) > 2.5)
                if operation == 'charge':
                    atoms['charge'][index] = str(float(atoms['charge'][index]) + .1)
                else:
                    atoms['chemical_type'][index] = ''
                document[0].set_mmcif_category('_cns_mtf_atom.', atoms)
            path.write_text(document.as_string())
        return directory, target

    def test_lost_or_moved_heavy_atom_is_rejected(self):
        for mutation, expected in [('missing_atom', 'complete_named_heavy_atom_inventory'),
                                   ('moved_atom', 'heavy_atom_coordinates_changed_or_rebuilt')]:
            directory, path = self._mutated_case(mutation)
            with directory:
                self.assertIn(expected, audit_topology('2RDK', path, root=self.root)['failures'])

    def test_actual_link_loss_charge_or_untyped_atom_is_rejected(self):
        for mutation, expected in [('missing_bond', 'complete_heavy_bond_connectivity'),
                                   ('charge', 'parent_total_charge'),
                                   ('type', 'missing_CNS_heavy_atom_type')]:
            directory, path = self._mutated_case(mutation)
            with directory:
                self.assertIn(expected, audit_topology('2RDK', path, root=self.root)['failures'])

    def test_source_byte_tampering_fails_before_conversion(self):
        with tempfile.TemporaryDirectory(dir=self.output) as directory:
            path = Path(directory) / 'source'
            path.write_bytes(b'original')
            expected = hashlib.sha256(path.read_bytes()).hexdigest()
            path.write_bytes(b'tampered')
            with self.assertRaisesRegex(CohortPreparationError, 'checksum'):
                _checked(path, expected)


if __name__ == '__main__':
    unittest.main()
