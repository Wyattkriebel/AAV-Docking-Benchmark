#!/usr/bin/env python3
"""Source-only strict/repair HADDOCK topology preflight; never docks."""
from collections import Counter
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
from rdkit import Chem

from aav_benchmark.haddock_cohort_preparation import prepare_glycan, audit_topology, pin
from aav_benchmark.haddock_method import pdb_atoms, read_cns_topology

POLICY = ROOT / 'benchmarks/haddock_paired_topology_policy.json'
OUTPUT = ROOT / 'results/haddock_prospective_preparation/paired_receptors'


def _key(atom):
    return tuple(atom[k] for k in ['chain', 'resnum', 'insertion', 'resname', 'atom_name'])


def prepare_pair(case, profile, output=OUTPUT):
    source = ROOT / 'data/processed/heldout_pose_v3' / profile / case
    destination = output / profile / case
    if (destination / 'preparation.json').exists():
        raise ValueError('Refuse to overwrite paired preflight')
    prepared = prepare_glycan(case, destination)
    receptor = source / 'receptor.pdb'
    protein = pdb_atoms(receptor)
    if 'L' in {a['chain'] for a in protein}:
        raise ValueError('Reserved ligand chain L collides with receptor')
    shutil.copyfile(receptor, destination / 'data/receptor.pdb')
    ligand = pdb_atoms(destination / 'data/glycan.pdb')
    translation = np.array([max(a['xyz_A'][0] for a in protein) + 50 - min(a['xyz_A'][0] for a in ligand), 0, 0])
    lines = (destination / 'data/glycan.pdb').read_text().splitlines()
    converted = []
    for line in lines:
        if line.startswith(('ATOM  ', 'HETATM')):
            xyz = np.array([float(line[start:start + 8]) for start in [30, 38, 46]]) + translation
            line = line[:21] + 'L' + line[22:30] + ''.join(f'{x:8.3f}' for x in xyz) + line[54:72] + 'L   ' + line[76:]
        converted.append(line)
    (destination / 'data/glycan.pdb').write_text('\n'.join(converted) + '\n')
    sdf_source = ROOT / prepared['input_sdf']['path']
    molecule = Chem.SDMolSupplier(str(sdf_source), removeHs=False)[0]
    conformer = molecule.GetConformer()
    for i in range(molecule.GetNumAtoms()):
        conformer.SetAtomPosition(i, np.array(conformer.GetAtomPosition(i)) + translation)
    sdf = destination / 'data/glycan_translated.sdf'
    writer = Chem.SDWriter(str(sdf)); writer.write(molecule); writer.close()
    prepared.update(original_independent_sdf=prepared['input_sdf'], input_sdf=pin(sdf),
                    independent_rigid_translation_A=translation.tolist(), receptor_source=pin(receptor),
                    observed_protein_atom_map=pin(source / 'protein_atom_map.json'))
    if (source / 'modeled_atom_map.json').exists():
        prepared['modeled_repair_atom_map_not_truth'] = pin(source / 'modeled_atom_map.json')
    for atom in prepared['atom_map']:
        atom['cns_chain'] = 'L'
    prepared['input_PDB'] = pin(destination / 'data/glycan.pdb')
    config = ('run_dir = "run-topology"\nmode = "local"\nncores = 1\ndebug = true\nclean = false\n'
              'molecules = ["data/receptor.pdb", "data/glycan.pdb"]\n[topoaa]\n')
    (destination / 'topology.cfg').write_text(config)
    prepared['config'] = pin(destination / 'topology.cfg')
    (destination / 'preparation.json').write_text(json.dumps(prepared, indent=2) + '\n')
    return destination, prepared


def audit_receptor(destination, policy):
    input_pdb = destination / 'data/receptor.pdb'
    output_pdb = destination / 'run-topology/0_topoaa/receptor_haddock.pdb'
    topology = destination / 'run-topology/0_topoaa/receptor_haddock.psf'
    if not output_pdb.is_file() or not topology.is_file():
        return dict(status='failed_no_receptor_topology')
    before = {_key(a): a for a in pdb_atoms(input_pdb) if a['element'] not in ['H', 'D']}
    after = {_key(a): a for a in pdb_atoms(output_pdb) if a['element'] not in ['H', 'D']}
    shared = sorted(set(before) & set(after))
    changes = [dict(atom=list(k), displacement_A=float(np.linalg.norm(np.array(before[k]['xyz_A']) - after[k]['xyz_A']))) for k in shared]
    changed = [a for a in changes if a['displacement_A'] > policy['protein_coordinate_tolerance_A']]
    topo = read_cns_topology(topology)
    residues = {}
    for key in shared:
        if key[4] in ['N', 'CA', 'C', 'CB']:
            residues.setdefault(key[:4], {})[key[4]] = key
    stereo = []
    for residue, atoms in residues.items():
        if set(atoms) != {'N', 'CA', 'C', 'CB'}:
            continue
        def volume(model):
            center = np.array(model[atoms['CA']]['xyz_A'])
            return float(np.linalg.det(np.array([np.array(model[atoms[n]]['xyz_A']) - center for n in ['N', 'C', 'CB']])))
        v0, v1 = volume(before), volume(after)
        stereo.append(dict(residue=list(residue), input_volume=v0, output_volume=v1,
                           preserved=v0 * v1 > 0 and abs(v1) >= policy['protein_chiral_signed_volume_min_A3']))
    exact = set(before) == set(after) and not changed and all(s['preserved'] for s in stereo)
    return dict(status='exact_named_heavy_receptor_preserved' if exact else 'CNS_receptor_preparation_differs',
                input_heavy_atoms=len(before), output_heavy_atoms=len(after), shared_heavy_atoms=len(shared),
                input_chains=sorted({k[0] for k in before}), output_chains=sorted({k[0] for k in after}),
                input_only_heavy_atoms=[list(k) for k in sorted(set(before) - set(after))],
                CNS_added_or_renamed_heavy_atoms=[list(k) for k in sorted(set(after) - set(before))],
                shared_heavy_max_change_A=max((c['displacement_A'] for c in changes), default=None), changed_shared_heavy_atoms=changed,
                CA_stereo=stereo, CNS_total_partial_charge=topo['total_charge'], CNS_atom_types=dict(Counter(a['atom_type'] for a in topo['atoms'].values())),
                protein_charge_interpretation='Recorded actual CNS protonation/partial charge; no unsupported exact total formal-charge equivalence is asserted.',
                prepared_pdb=pin(output_pdb), actual_CNS_topology=pin(topology))


def main():
    policy = json.loads(POLICY.read_text())
    jobs = []
    for profile in policy['receptor_profiles']:
        for case in policy['cases']:
            destination, prepared = prepare_pair(case, profile)
            jobs.append((profile, case, destination, prepared))
    plan = dict(schema_version=1, source_only=True, no_docking_stages=True, written_before_first_CNS_job=True,
                policy=pin(POLICY), driver=pin(Path(__file__)),
                jobs=[dict(profile=p, case_id=c, receptor=prep['receptor_source'], independent_sdf=prep['original_independent_sdf'],
                           translation=prep['independent_rigid_translation_A'], input_PDB=prep['input_PDB'], config=prep['config'],
                           preparation=pin(d / 'preparation.json')) for p, c, d, prep in jobs])
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / 'execution_plan.json').write_text(json.dumps(plan, indent=2) + '\n')
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1', MKL_NUM_THREADS='1')
    records = []
    for profile, case, destination, prepared in jobs:
        start = time.monotonic()
        with (destination / 'stdout.log').open('w') as stdout, (destination / 'stderr.log').open('w') as stderr:
            result = subprocess.run([str(ROOT / 'tools/haddock3-2026.10.1/env/bin/haddock3'), 'topology.cfg'], cwd=destination,
                                    env=environment, stdout=stdout, stderr=stderr, check=False)
        record = dict(profile=profile, case_id=case, source_only=True, no_docking_performed=True,
                      process_return_code=result.returncode, elapsed_seconds=time.monotonic() - start,
                      glycan=audit_topology(case, destination), receptor=audit_receptor(destination, policy))
        (destination / 'paired_topology_audit.json').write_text(json.dumps(record, indent=2) + '\n')
        records.append(record)
        print(profile, case, result.returncode, record['glycan']['status'], record['receptor']['status'], flush=True)
    (OUTPUT / 'summary.json').write_text(json.dumps(dict(schema_version=1, source_only=True, no_docking_performed=True,
                                                       policy=pin(POLICY), driver=pin(Path(__file__)), cases=records), indent=2) + '\n')


if __name__ == '__main__':
    main()
