"""Chemical-only controls for the separate ideal-ring initialization profile."""
from copy import deepcopy
import gzip
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import gemmi
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

from aav_benchmark import ideal_ring_conformers as ir
from aav_benchmark.ring_geometry import analytic_ring, compare_to_ideal, cremer_pople6


class IdealRingControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.policy = json.loads((ir.ROOT / ir.POLICY_PATH).read_text())
        cls.parents = {case: ir.build_parent_graph(ir.case_schema(case)) for case in cls.policy['cases']}
        cls.output = ir.ROOT / 'results/ideal_ring_conformer_feasibility' / cls.policy['profile_id']
        cls.selected = Chem.SDMolSupplier(str(cls.output / '1GNY/ligand.sdf'), removeHs=False)[0]

    def test_all_seven_complete_chemical_graphs_equal_independent_frozen_graphs(self):
        golden = json.loads((ir.ROOT / 'tests/fixtures/ideal_ring_expected_parent_graphs.json').read_text())
        self.assertEqual(set(golden['cases']), set(self.parents))
        for case, (mol, mapping, rings, stereo, provenance) in self.parents.items():
            with self.subTest(case=case):
                self.assertEqual(provenance['isomeric_smiles'], golden['cases'][case]['canonical_isomeric_smiles'])
                self.assertEqual(mol.GetNumConformers(), 0)
                self.assertFalse(provenance['experimental_coordinates_read'])
                self.assertEqual(len(mapping), mol.GetNumHeavyAtoms())
                self.assertEqual(len(set(a['canonical_atom_key'] for a in mapping)), len(mapping))
                self.assertEqual(len(Chem.GetMolFrags(mol)), 1)
                self.assertEqual(ir.graph_fingerprint(mol), mol.GetProp('expected_chemical_graph'))
                self.assertEqual(len(stereo), len(Chem.FindMolChiralCenters(Chem.RemoveHs(mol), includeUnassigned=True)))

    def test_builder_reads_only_pinned_chemical_sources(self):
        schema = ir.case_schema('1GNY')
        allowed = {(ir.ROOT / (s.get('archive') or s['archive_path'])).resolve() for s in schema['ccd_sources'].values()}
        original = Path.read_bytes
        observed = []
        def guarded(path):
            self.assertIn(path.resolve(), allowed)
            observed.append(path.resolve())
            return original(path)
        with patch.object(Path, 'read_bytes', guarded):
            ir.build_parent_graph(schema)
        self.assertEqual(set(observed), allowed)

    def test_each_chemical_ideal_component_agrees_with_declared_ccd_stereochemistry(self):
        seen = set()
        for case in self.policy['cases']:
            for name, source in ir.case_schema(case)['ccd_sources'].items():
                if name in seen:
                    continue
                seen.add(name)
                component = ir._component(source, ir.ROOT)
                for atom_name, definition in component['atoms'].items():
                    if definition['stereo'] in ('R', 'S'):
                        atom = component['mol'].GetAtomWithIdx(component['indices'][atom_name])
                        self.assertEqual(atom.GetProp('_CIPCode'), definition['stereo'])

    def test_experimental_ccd_model_coordinates_are_ignored(self):
        schema = ir.case_schema('1GNY')
        source = next(iter(schema['ccd_sources'].values()))
        original = ir._component(source, ir.ROOT)
        path = ir.ROOT / (source.get('archive') or source['archive_path'])
        raw = path.read_bytes()
        payload = gzip.decompress(raw) if path.suffix == '.gz' else raw
        document = gemmi.cif.read_string(payload.decode())
        atoms = document[0].get_mmcif_category('_chem_comp_atom.')
        for axis in 'xyz':
            atoms['model_Cartn_' + axis] = [str(99999 + i * 1234) for i in range(len(atoms['atom_id']))]
        document[0].set_mmcif_category('_chem_comp_atom.', atoms)
        poisoned = document.as_string().encode()
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / 'poisoned.cif.gz'
            target.write_bytes(gzip.compress(poisoned, mtime=0))
            replacement = deepcopy(source)
            replacement.update(archive=str(target), archive_path=str(target),
                               archive_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
                               sha256=hashlib.sha256(poisoned).hexdigest())
            actual = ir._component(replacement, ir.ROOT)
        np.testing.assert_array_equal(actual['ring_xyz'], original['ring_xyz'])
        np.testing.assert_array_equal(actual['ring_torsions'], original['ring_torsions'])
        self.assertEqual(Chem.MolToSmiles(actual['mol']), Chem.MolToSmiles(original['mol']))

    def _validation(self, mol):
        _, mapping, rings, stereo, _ = self.parents['1GNY']
        return ir.validate_conformer(mol, 0, mapping, rings, stereo, self.policy)

    def test_selected_serialized_sdf_passes_full_chemical_and_physical_gates(self):
        self.assertTrue(self._validation(self.selected)['valid'])

    def test_bond_charge_chiral_tag_and_atom_loss_cannot_pass_graph_gate(self):
        variants = []
        mol = Chem.RWMol(self.selected)
        mol.RemoveBond(0, 1)
        variants.append(mol.GetMol())
        mol = Chem.RWMol(self.selected)
        mol.GetBondWithIdx(0).SetBondType(Chem.BondType.DOUBLE)
        variants.append(mol.GetMol())
        mol = Chem.Mol(self.selected)
        mol.GetAtomWithIdx(0).SetFormalCharge(1)
        variants.append(mol)
        mol = Chem.Mol(self.selected)
        next(a for a in mol.GetAtoms() if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED).InvertChirality()
        variants.append(mol)
        mol = Chem.RWMol(self.selected)
        mol.RemoveAtom(0)
        variants.append(mol.GetMol())
        for mol in variants:
            try:
                result = self._validation(mol)
            except (ValueError, RuntimeError, KeyError):
                continue
            self.assertFalse(result['valid'])
            self.assertIn('full_parent_bond_charge_or_isomeric_graph_mismatch', result['failures'])

    def test_reflected_and_planar_named_stereocenters_fail(self):
        _, mapping, _, stereo, _ = self.parents['1GNY']
        keys = {a['canonical_atom_key']: a['rdkit_index'] for a in mapping}
        center = next(s for s in stereo if len(s['neighbors']) == 3)
        neighbors = self.selected.GetConformer().GetPositions()[[keys[n] for n in center['neighbors']]]
        normal = np.cross(neighbors[1] - neighbors[0], neighbors[2] - neighbors[0])
        normal /= np.linalg.norm(normal)
        initial = np.array(self.selected.GetConformer().GetAtomPosition(keys[center['center']]))
        reflected = initial - 2 * np.dot(initial - neighbors[0], normal) * normal
        planar = neighbors.mean(axis=0)
        for position in (reflected, planar):
            mol = Chem.Mol(self.selected)
            mol.GetConformer().SetAtomPosition(keys[center['center']], position)
            self.assertIn('full_parent_stereochemistry', self._validation(mol)['failures'])

    def test_opposite_chair_and_boat_are_rejected_by_named_pucker(self):
        south, north, boat = analytic_ring(0, -.6), analytic_ring(0, .6), analytic_ring(.6, 0)
        self.assertLess(cremer_pople6(north)['chair_pole_distance_deg'], 1e-9)
        self.assertGreater(compare_to_ideal(north, south)['pucker_vector_angle_to_CCD_ideal_deg'], 179)
        self.assertGreater(compare_to_ideal(boat, south)['pucker_vector_angle_to_CCD_ideal_deg'], 89)

    def test_fixed_slot_budget_coordinate_discard_and_exclusively_intraring_restraints(self):
        mol, mapping, rings, stereo, _ = self.parents['1GNY']
        poisoned = Chem.Mol(mol)
        conf = Chem.Conformer(poisoned.GetNumAtoms())
        for i in range(poisoned.GetNumAtoms()):
            conf.SetAtomPosition(i, (i * 10000, 50000, -50000))
        poisoned.AddConformer(conf)
        seeds, fields = [], []
        class Field:
            def __init__(self): self.distances, self.torsions = [], []
            def MMFFAddDistanceConstraint(self, *args): self.distances.append(args)
            def MMFFAddTorsionConstraint(self, *args): self.torsions.append(args)
            def Initialize(self): pass
            def Minimize(self, **kwargs): return 0
            def CalcEnergy(self): return 1.0
        def embed(candidate, params):
            self.assertEqual(candidate.GetNumConformers(), 0)
            self.assertEqual(params.numThreads, 1)
            seeds.append(params.randomSeed)
            candidate.AddConformer(Chem.Conformer(self.selected.GetConformer()))
            return 0
        def forcefield(*args, **kwargs):
            field = Field(); fields.append(field); return field
        with patch.object(ir.AllChem, 'EmbedMolecule', embed), patch.object(ir.AllChem, 'MMFFGetMoleculeForceField', forcefield):
            selected, result = ir.generate_conformer(poisoned, mapping, rings, stereo, self.policy)
        self.assertEqual(seeds, list(range(20261024, 20261044)))
        self.assertEqual(len(result['attempts']), 20)
        self.assertEqual(result['selected_slot'], 0)
        np.testing.assert_array_equal(selected.GetConformer().GetPositions(), self.selected.GetConformer().GetPositions())
        ring_sets = [set(r['indices']) for r in rings]
        constrained = [f for f in fields if f.distances]
        self.assertEqual(len(constrained), 20)
        for field in constrained:
            self.assertEqual(len(field.distances), len(rings) * 15)
            self.assertEqual(len(field.torsions), len(rings) * 6)
            for constraint in field.distances:
                self.assertTrue(any(set(constraint[:2]) <= r for r in ring_sets))
                self.assertFalse(constraint[2])
            for constraint in field.torsions:
                self.assertTrue(any(set(constraint[:4]) <= r for r in ring_sets))
                self.assertFalse(constraint[4])

    def test_glycosidic_rotation_preserves_all_intraring_restraint_observables(self):
        mol = Chem.Mol(self.selected)
        _, mapping, rings, _, provenance = self.parents['1GNY']
        keys = {a['canonical_atom_key']: a['rdkit_index'] for a in mapping}
        left, right = [keys[k] for k in provenance['declared_links'][0]]
        bond = mol.GetBondBetweenAtoms(left, right)
        disconnected = Chem.FragmentOnBonds(mol, [bond.GetIdx()], addDummies=False)
        moving = next(set(f) for f in Chem.GetMolFrags(disconnected) if left in f)
        before = mol.GetConformer().GetPositions().copy()
        origin, axis = before[right], before[left] - before[right]
        axis /= np.linalg.norm(axis)
        angle = .7
        for index in moving:
            p = before[index] - origin
            q = p * np.cos(angle) + np.cross(axis, p) * np.sin(angle) + axis * np.dot(axis, p) * (1 - np.cos(angle))
            mol.GetConformer().SetAtomPosition(index, q + origin)
        after = mol.GetConformer().GetPositions()
        for ring in rings:
            indices = ring['indices']
            np.testing.assert_allclose(np.linalg.norm(before[indices, None] - before[indices], axis=2),
                                       np.linalg.norm(after[indices, None] - after[indices], axis=2), atol=1e-12)
            for k in range(6):
                atoms = [indices[(k + q) % 6] for q in range(4)]
                self.assertLess(ir._angle_delta(rdMolTransforms.GetDihedralDeg(self.selected.GetConformer(), *atoms),
                                               rdMolTransforms.GetDihedralDeg(mol.GetConformer(), *atoms)), 1e-9)
        a = next(n.GetIdx() for n in mol.GetAtomWithIdx(left).GetNeighbors() if n.GetIdx() != right)
        b = next(n.GetIdx() for n in mol.GetAtomWithIdx(right).GetNeighbors() if n.GetIdx() != left)
        self.assertGreater(ir._angle_delta(rdMolTransforms.GetDihedralDeg(self.selected.GetConformer(), a, left, right, b),
                                          rdMolTransforms.GetDihedralDeg(mol.GetConformer(), a, left, right, b)), 35)

    def test_all_seven_saved_products_revalidate_and_selection_is_minimum_valid_energy(self):
        for case in self.policy['cases']:
            with self.subTest(case=case):
                report = json.loads((self.output / case / 'report.json').read_text())
                mol = Chem.SDMolSupplier(str(self.output / case / 'ligand.sdf'), removeHs=False)[0]
                _, mapping, rings, stereo, _ = self.parents[case]
                result = ir.validate_conformer(mol, 0, mapping, rings, stereo, self.policy)
                self.assertTrue(result['valid'], result['failures'])
                from_3d = Chem.Mol(mol)
                Chem.AssignAtomChiralTagsFromStructure(from_3d, replaceExistingTags=True)
                Chem.AssignStereochemistry(from_3d, cleanIt=True, force=True)
                self.assertEqual(Chem.MolToSmiles(Chem.RemoveHs(from_3d), isomericSmiles=True),
                                 self.parents[case][4]['isomeric_smiles'])
                valid = [a for a in report['search']['attempts'] if a['status'] == 'valid']
                best = min(valid, key=lambda a: (a['unconstrained_mmff_energy_at_selected_geometry'], a['slot']))
                self.assertEqual(report['search']['selected_slot'], best['slot'])
                self.assertEqual(len(report['search']['attempts']), 20)
                self.assertEqual(report['policy_sha256'], ir._sha(ir.ROOT / ir.POLICY_PATH))
                self.assertEqual(report['algorithm_source_sha256'], ir._sha(ir.ROOT / 'aav_benchmark/ideal_ring_conformers.py'))
                self.assertEqual(report['selected_sdf_sha256'], ir._sha(self.output / case / 'ligand.sdf'))
                self.assertTrue(report['no_docking_or_accuracy_evaluation'])


if __name__ == '__main__':
    unittest.main()
