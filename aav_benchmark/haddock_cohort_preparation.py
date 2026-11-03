"""New source-only HADDOCK cohort parameterization gates; no docking stages.

The released-example lock and helper are read without modification. Inputs are
the frozen independent chemical-ideal v3 full parents, never bound glycans.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import os

import numpy as np
from rdkit import Chem

from . import ideal_ring_conformers as ideal
from .haddock_method import pdb_atoms, read_cns_topology

ROOT = Path(__file__).resolve().parents[1]
POLICY = Path('benchmarks/haddock_prospective_policy.json')
DEFAULT_OUTPUT = Path('results/haddock_prospective_preparation')


class CohortPreparationError(ValueError):
    pass


def pin(path, root=ROOT):
    path = Path(path)
    return dict(path=str(path.relative_to(root)), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                bytes=path.stat().st_size)


def _checked(path, expected):
    if hashlib.sha256(Path(path).read_bytes()).hexdigest() != expected:
        raise CohortPreparationError('Source checksum mismatch: ' + str(path))


def prepare_glycan(case_id, destination, root=ROOT):
    """Write standard CCD heavy names from v3, preserving its coordinates."""
    root, destination = Path(root), Path(destination)
    policy = json.loads((root / POLICY).read_text())
    if case_id not in policy['cases']:
        raise CohortPreparationError('Case outside fixed six-parent inventory')
    _checked(root / policy['independent_ligand_receipt'], policy['independent_ligand_receipt_sha256'])
    receipt = json.loads((root / policy['independent_ligand_receipt']).read_text())
    frozen = next(c for c in receipt['cases'] if c['case_id'] == case_id)
    sdf = root / frozen['ligand_sdf']['path']
    report_path = root / frozen['report']['path']
    _checked(sdf, frozen['ligand_sdf']['sha256'])
    _checked(report_path, frozen['report']['sha256'])
    report = json.loads(report_path.read_text())
    molecule = Chem.SDMolSupplier(str(sdf), removeHs=False)[0]
    graph, mapping, rings, stereo, provenance = ideal.build_parent_graph(ideal.case_schema(case_id, root), root)
    validation = ideal.validate_conformer(molecule, 0, mapping, rings, stereo,
                                         json.loads((root / ideal.POLICY_PATH).read_text()))
    if not validation['valid']:
        raise CohortPreparationError('Independent full parent fails physical/chemical gate')
    lines, canonical_map = [], []
    for serial, atom in enumerate(mapping, 1):
        xyz = molecule.GetConformer().GetAtomPosition(atom['rdkit_index'])
        number, name, comp = int(atom['branch_num']), atom['atom_name'], atom['ccd_id']
        if len(name) > 4 or len(comp) > 3 or number > 9999:
            raise CohortPreparationError('PDB/CNS name or residue-number overflow')
        line = (f"HETATM{serial:5d} {name:>4s} {comp:3s} B{number:4d}    "
                f"{xyz.x:8.3f}{xyz.y:8.3f}{xyz.z:8.3f}{1.0:6.2f}{0.0:6.2f}      B   {atom['element']:>2s}  ")
        lines.append(line)
        canonical_map.append(dict(**atom, cns_chain='B', cns_resnum=str(number), cns_resname=comp,
                                  cns_atom_name=name, input_pdb_serial=serial))
    destination.mkdir(parents=True, exist_ok=True)
    data = destination / 'data'
    data.mkdir(exist_ok=True)
    pdb = data / 'glycan.pdb'
    pdb.write_text('\n'.join(lines) + '\nTER\nEND\n')
    if len(pdb_atoms(pdb)) != len(mapping):
        raise CohortPreparationError('PDB conversion lost atoms')
    config = ('run_dir = "run-topology"\nmode = "local"\nncores = 1\n'
              'debug = true\nclean = false\nmolecules = ["data/glycan.pdb"]\n[topoaa]\n')
    (destination / 'topology.cfg').write_text(config)
    prepared = dict(schema_version=1, case_id=case_id, source_only=True, native_coordinates_read=False,
                    no_docking_stages=True, input_sdf=pin(sdf, root), source_report=pin(report_path, root),
                    full_parent_chemistry=provenance, atom_map=canonical_map,
                    expected_heavy_bonds=[dict(left=graph.GetAtomWithIdx(b.GetBeginAtomIdx()).GetProp('canonical_atom_key'),
                                              right=graph.GetAtomWithIdx(b.GetEndAtomIdx()).GetProp('canonical_atom_key'),
                                              order=str(b.GetBondType()))
                                          for b in graph.GetBonds()
                                          if b.GetBeginAtom().GetAtomicNum() != 1 and b.GetEndAtom().GetAtomicNum() != 1],
                    input_PDB=pin(pdb, root), config=pin(destination / 'topology.cfg', root),
                    hydrogen_policy='Input heavy atoms only; all output hydrogens are explicitly CNS-built.')
    (destination / 'preparation.json').write_text(json.dumps(prepared, indent=2) + '\n')
    return prepared


def audit_topology(case_id, destination, root=ROOT):
    """Require actual complete CNS topology and geometry, not topology names."""
    root, destination = Path(root), Path(destination)
    policy = json.loads((root / POLICY).read_text())
    prepared = json.loads((destination / 'preparation.json').read_text())
    stage = destination / 'run-topology/0_topoaa'
    pdb = stage / 'glycan_haddock.pdb'
    topology = stage / 'glycan_haddock.psf'
    if not pdb.is_file() or not topology.is_file():
        return dict(case_id=case_id, status='failed_no_prepared_topology',
                    stage_files=[p.name for p in stage.glob('*')] if stage.is_dir() else [])
    before = {(a['resnum'], a['resname'], a['atom_name']): a for a in pdb_atoms(destination / 'data/glycan.pdb')}
    after_all = pdb_atoms(pdb)
    after = {(a['resnum'], a['resname'], a['atom_name']): a for a in after_all if a['element'] not in ['H', 'D']}
    declared = {(a['cns_resnum'], a['cns_resname'], a['cns_atom_name']): a['canonical_atom_key'] for a in prepared['atom_map']}
    topo = read_cns_topology(topology)
    topo_heavy = {i: a for i, a in topo['atoms'].items() if a['mass'] > 2.5}
    topo_keys = {(a['resnum'], a['resname'], a['atom_name']) for a in topo_heavy.values()}
    failures = []
    if set(before) != set(after) or topo_keys != set(declared):
        failures.append('complete_named_heavy_atom_inventory')
    actual_edges = set()
    for left, right in topo['bonds']:
        if left in topo_heavy and right in topo_heavy:
            keys = [tuple(topo_heavy[i][k] for k in ['resnum', 'resname', 'atom_name']) for i in [left, right]]
            if all(k in declared for k in keys):
                actual_edges.add(tuple(sorted(declared[k] for k in keys)))
            else:
                failures.append('unexpected_heavy_atom_in_bond')
    expected_edges = {tuple(sorted((e['left'], e['right']))) for e in prepared['expected_heavy_bonds']}
    if actual_edges != expected_edges:
        failures.append('complete_heavy_bond_connectivity')
    if any(e['order'] != 'SINGLE' for e in prepared['expected_heavy_bonds']):
        failures.append('bond_order_audit_not_implemented_for_non_single_parent')
    if any(not a['atom_type'] or a['mass'] <= 0 for a in topo_heavy.values()):
        failures.append('missing_CNS_heavy_atom_type')
    charge_delta = topo['total_charge'] - prepared['full_parent_chemistry']['formal_charge']
    if abs(charge_delta) > policy['topology_acceptance']['total_CNS_partial_charge_matches_parent_formal_charge_tolerance_e']:
        failures.append('parent_total_charge')
    displacements = {':'.join(k): float(np.linalg.norm(np.array(before[k]['xyz_A']) - after[k]['xyz_A']))
                     for k in set(before) & set(after)}
    if max(displacements.values(), default=float('inf')) > policy['topology_acceptance']['max_named_heavy_coordinate_change_from_input_PDB_A']:
        failures.append('heavy_atom_coordinates_changed_or_rebuilt')
    geometry, from_3d_smiles = None, None
    if set(declared) == set(after):
        molecule = Chem.SDMolSupplier(str(root / prepared['input_sdf']['path']), removeHs=False)[0]
        _, mapping, rings, stereo, _ = ideal.build_parent_graph(ideal.case_schema(case_id, root), root)
        indices = {a['canonical_atom_key']: a['rdkit_index'] for a in mapping}
        for key, canonical in declared.items():
            molecule.GetConformer().SetAtomPosition(indices[canonical], after[key]['xyz_A'])
        geometry = ideal.validate_conformer(molecule, 0, mapping, rings, stereo,
                                            json.loads((root / ideal.POLICY_PATH).read_text()))
        if not geometry['valid']:
            failures.append('post_CNS_full_parent_geometry_or_stereo')
        Chem.AssignAtomChiralTagsFromStructure(molecule, replaceExistingTags=True)
        Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
        from_3d_smiles = Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True)
        if from_3d_smiles != prepared['full_parent_chemistry']['isomeric_smiles']:
            failures.append('post_CNS_3D_isomeric_graph')
    return dict(schema_version=1, case_id=case_id, status='accepted_glycan_parameterization' if not failures else 'rejected_glycan_parameterization',
                source_only=True, no_docking_performed=True, failures=sorted(set(failures)),
                prepared_pdb=pin(pdb, root), actual_CNS_topology=pin(topology, root),
                input_heavy_count=len(before), output_heavy_count=len(after), typed_heavy_count=len(topo_heavy),
                CNS_hydrogens_built=sum(a['mass'] < 2.5 for a in topo['atoms'].values()),
                topology_total_charge=topo['total_charge'], expected_parent_formal_charge=prepared['full_parent_chemistry']['formal_charge'],
                total_charge_difference=charge_delta, heavy_coordinate_max_change_A=max(displacements.values(), default=None),
                missing_input_heavy_atoms=[list(k) for k in sorted(set(before) - set(after))],
                extra_output_heavy_atoms=[list(k) for k in sorted(set(after) - set(before))],
                missing_heavy_bonds=[list(k) for k in sorted(expected_edges - actual_edges)],
                extra_heavy_bonds=[list(k) for k in sorted(actual_edges - expected_edges)],
                actual_heavy_bonds=[list(k) for k in sorted(actual_edges)],
                post_CNS_geometry=geometry, fresh_3D_isomeric_smiles=from_3d_smiles,
                interpretation='Actual glycan-only preparation gate. Protein/receptor topology and box-derived AIR fairness are separate gates before any cohort docking.')


def run_parameterization(root=ROOT, output=None, cases=None, execute=False):
    root = Path(root)
    output = Path(output or root / DEFAULT_OUTPUT)
    policy = json.loads((root / POLICY).read_text())
    lock = json.loads((root / policy['release_lock_path']).read_text())
    for key in ['distribution', 'cns_binary']:
        _checked(root / lock[key]['path'], lock[key]['sha256'])
    selected_cases = cases or policy['cases']
    input_plans = []
    for case_id in selected_cases:
        destination = output / case_id
        if (destination / 'preparation.json').exists():
            raise CohortPreparationError('Refuse to overwrite an existing preparation')
        preparation = prepare_glycan(case_id, destination, root)
        input_plans.append(dict(case_id=case_id, input_PDB=preparation['input_PDB'],
                                config=preparation['config'], preparation=pin(destination / 'preparation.json', root)))
    plan = dict(schema_version=1, source_only=True, no_docking_stages=True,
                written_before_first_CNS_job=True, execute_topoaa_only=execute,
                policy=pin(root / POLICY, root), algorithm=pin(Path(__file__), root),
                CNS_binary=lock['cns_binary'], inputs=input_plans)
    output.mkdir(parents=True, exist_ok=True)
    (output / 'execution_plan.json').write_text(json.dumps(plan, indent=2) + '\n')
    records = []
    for case_id in selected_cases:
        destination = output / case_id
        if execute:
            environment = os.environ.copy()
            environment.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1', MKL_NUM_THREADS='1')
            command = [str(root / 'tools/haddock3-2026.10.1/env/bin/haddock3'), 'topology.cfg']
            with (destination / 'stdout.log').open('w') as stdout, (destination / 'stderr.log').open('w') as stderr:
                completed = subprocess.run(command, cwd=destination, env=environment, stdout=stdout, stderr=stderr, check=False)
            record = audit_topology(case_id, destination, root)
            record['process_return_code'] = completed.returncode
            if completed.returncode:
                record['status'] = 'failed_CNS_preparation'
            (destination / 'topology_audit.json').write_text(json.dumps(record, indent=2) + '\n')
        else:
            record = dict(case_id=case_id, status='prepared_not_executed')
        records.append(record)
        print(case_id, record['status'], flush=True)
    summary = dict(schema_version=1, source_only=True, no_docking_performed=True,
                   policy=pin(root / POLICY, root), algorithm=pin(Path(__file__), root),
                   CNS_binary=lock['cns_binary'], cases=records)
    output.mkdir(parents=True, exist_ok=True)
    (output / 'parameterization_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--cases', nargs='+')
    parser.add_argument('--execute-topoaa-only', action='store_true')
    args = parser.parse_args(argv)
    summary = run_parameterization(args.root, args.output, args.cases, args.execute_topoaa_only)
    return 0 if all(not c.get('failures') and not c['status'].startswith('failed') for c in summary['cases']) else 1


if __name__ == '__main__':
    raise SystemExit(main())
