"""Audit a pinned, actually executed HADDOCK3 released-example reproduction.

This is a known-binding-site guided software reproduction, not a blind benchmark.
No source qualification ledger is changed and no new cohort case is executed.
"""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import re
import tarfile

import gemmi
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_LOCK = Path('benchmarks/haddock.lock.json')
DEFAULT_EXAMPLE = Path('results/haddock_reproduction/official_glycan_test')
EXPECTED_SDIST_SHA256 = '2edf1f7efef5d9d76217f0c378c14781461a5f729c56ad4b91b1926f0d15421c'
EXPECTED_CNS_SHA256 = 'e6febf74eccb740f8bc5d49955cc8fc202025d5d642e129e7738e7280ad95a88'
GLYCAM_ACETYL_ALIASES = {'C2N': 'C7', 'O2N': 'O7', 'CME': 'C8'}


class HaddockAuditError(ValueError):
    """A release, chemical inventory or completed workflow invariant failed."""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _text(path: Path) -> str:
    return gzip.open(path, 'rt').read() if path.suffix == '.gz' else path.read_text()


def pdb_atoms(path: Path | str) -> list[dict]:
    """Read every PDB atom without inventing missing atoms or alias conversion."""
    atoms = []
    for line in _text(Path(path)).splitlines():
        if not line.startswith(('ATOM  ', 'HETATM')):
            continue
        name = line[12:16].strip()
        element = line[76:78].strip().upper()
        if not element:
            element = next((c.upper() for c in name if c.isalpha()), '')
        atoms.append(dict(chain=line[21:22].strip(), resnum=line[22:26].strip(),
                          insertion=line[26:27].strip(), resname=line[17:20].strip(),
                          atom_name=name, element=element, alt=line[16:17].strip(),
                          xyz_A=[float(line[start:start+8]) for start in [30, 38, 46]]))
    if not atoms:
        raise HaddockAuditError(f'No PDB atoms: {path}')
    keys = [(a['chain'], a['resnum'], a['insertion'], a['resname'], a['atom_name'], a['alt']) for a in atoms]
    if len(keys) != len(set(keys)):
        raise HaddockAuditError(f'Duplicate PDB atom identity: {path}')
    return atoms


def pdb_inventory(path: Path | str) -> dict:
    atoms = pdb_atoms(path)
    heavy = [a for a in atoms if a['element'] not in ['H', 'D']]
    residues = {}
    for atom in heavy:
        key = ':'.join(atom[k] for k in ['chain', 'resnum', 'resname'])
        residues.setdefault(key, []).append(atom['atom_name'])
    return dict(path=str(path), sha256=sha256(Path(path)), atom_count=len(atoms),
                heavy_atom_count=len(heavy), residue_count=len(residues),
                residue_heavy_atoms=residues, chain_heavy_atoms=dict(Counter(a['chain'] for a in heavy)))


def read_cns_topology(path: Path | str) -> dict:
    """Read the actual CNS mmCIF MTF written with the .psf extension."""
    block = gemmi.cif.read_string(_text(Path(path)))[0]
    category = block.get_mmcif_category('_cns_mtf_atom.')
    for i,index in enumerate(category['id']):
        if int(index) <= 0 and not (index == '-1' and
            all(not category[field][i].strip() for field in ['segment_id','residue_id','residue_name','atom_name','chemical_type']) and
            float(category['charge'][i]) == -1 and float(category['atom_mass'][i]) == -1):
            raise HaddockAuditError('Unexpected nonpositive CNS topology atom record')
    atoms = {category['id'][i]: dict(segment=category['segment_id'][i],
             resnum=category['residue_id'][i], resname=category['residue_name'][i],
             atom_name=category['atom_name'][i], atom_type=category['chemical_type'][i],
             charge=float(category['charge'][i]), mass=float(category['atom_mass'][i]))
             for i in range(len(category['id'])) if int(category['id'][i]) > 0}
    bonds = block.get_mmcif_category('_cns_mtf_bond.')
    edges = [(left, bonds['id[2]'][i]) for i, left in enumerate(bonds['id[1]'])
             if not (left == '-1' and bonds['id[2]'][i] == '-1')]
    if any(left not in atoms or right not in atoms for left, right in edges):
        raise HaddockAuditError('CNS topology bond references an absent atom')
    return dict(atoms=atoms, bonds=edges,
                total_charge=sum(a['charge'] for a in atoms.values()),
                atom_count=len(atoms), bond_count=len(edges),
                sentinel_atom_rows_ignored=len(category['id'])-len(atoms),
                sentinel_bond_rows_ignored=len(bonds['id[1]'])-len(edges))


def preparation_correspondence(input_pdb: Path | str, prepared_pdb: Path | str) -> dict:
    """Expose atom rebuilding; semantically mapped aliases are not called renames."""
    before = {(a['resnum'], a['resname'], a['atom_name']): a for a in pdb_atoms(input_pdb)
              if a['element'] not in ['H', 'D']}
    after = {(a['resnum'], a['resname'], a['atom_name']): a for a in pdb_atoms(prepared_pdb)
             if a['element'] not in ['H', 'D']}
    shared = sorted(set(before) & set(after))
    displacements = [float(np.linalg.norm(np.array(before[key]['xyz_A']) - after[key]['xyz_A'])) for key in shared]
    alias_comparisons = []
    for key, atom in before.items():
        if key[2] not in GLYCAM_ACETYL_ALIASES:
            continue
        output_key = (*key[:2], GLYCAM_ACETYL_ALIASES[key[2]])
        target = after.get(output_key)
        alias_comparisons.append(dict(input_atom=list(key), standard_semantic_atom=list(output_key),
            output_present=target is not None,
            coordinate_displacement_A=(float(np.linalg.norm(np.array(atom['xyz_A']) - target['xyz_A'])) if target else None)))
    return dict(input_heavy_atoms=len(before), prepared_heavy_atoms=len(after),
                exact_name_shared_heavy_atoms=len(shared),
                shared_max_displacement_A=max(displacements, default=0),
                input_only_atoms=[list(k) for k in sorted(set(before)-set(after))],
                prepared_only_atoms=[list(k) for k in sorted(set(after)-set(before))],
                semantic_alias_coordinate_comparisons=alias_comparisons,
                interpretation='Atom-name disappearance/addition and coordinates are observed. Semantic alias labels do not prove a coordinate-preserving rename; absent standard atoms may be rebuilt by CNS.')


def glycan_chirality(pdb: Path | str, topology: Path | str, ccd_sources: dict[str, Path]) -> dict:
    """Compare named CCD handedness with actual output coordinates and CNS bonds."""
    model = pdb_atoms(pdb)
    by_key = {(a['chain'], a['resnum'], a['resname'], a['atom_name']): a for a in model}
    actual_topology = read_cns_topology(topology)
    topo_atoms = actual_topology['atoms']
    by_topology_key = {(a['segment'], a['resnum'], a['resname'], a['atom_name']): index
                       for index, a in topo_atoms.items()}
    expected_heavy_keys={key for key,index in by_topology_key.items()
                         if key[2] in ccd_sources and topo_atoms[index]['mass'] > 2.5}
    observed_heavy_keys={key for key,a in by_key.items()
                         if key[2] in ccd_sources and a['element'] not in ['H','D']}
    if expected_heavy_keys != observed_heavy_keys:
        raise HaddockAuditError('Prepared glycan heavy-atom inventory differs from actual CNS topology')
    neighbors = {index: set() for index in topo_atoms}
    for left, right in actual_topology['bonds']:
        neighbors[left].add(right); neighbors[right].add(left)
    residue_keys = sorted({(a['chain'], a['resnum'], a['resname']) for a in model if a['resname'] in ccd_sources})
    centers = []
    links = []
    for left, right in actual_topology['bonds']:
        a, b = topo_atoms[left], topo_atoms[right]
        if (a['segment'], a['resnum']) != (b['segment'], b['resnum']):
            links.append(dict(left={k:a[k] for k in ['segment','resnum','resname','atom_name']},
                              right={k:b[k] for k in ['segment','resnum','resname','atom_name']}))
    for chain, resnum, resname in residue_keys:
        block = gemmi.cif.read_file(str(ccd_sources[resname]))[0]
        cc = block.get_mmcif_category('_chem_comp_atom.')
        bb = block.get_mmcif_category('_chem_comp_bond.')
        ideal = {name: np.array([float(cc['pdbx_model_Cartn_'+axis+'_ideal'][i]) for axis in 'xyz'])
                 for i,name in enumerate(cc['atom_id'])}
        elements = dict(zip(cc['atom_id'], cc['type_symbol']))
        edges = [(bb['atom_id_1'][i],bb['atom_id_2'][i]) for i in range(len(bb['atom_id_1']))
                 if elements[bb['atom_id_1'][i]] not in ['H','D'] and elements[bb['atom_id_2'][i]] not in ['H','D']]
        for i,name in enumerate(cc['atom_id']):
            if cc['pdbx_stereo_config'][i] not in ['R','S']:
                continue
            key=(chain,resnum,resname,name)
            if key not in by_key or key not in by_topology_key:
                raise HaddockAuditError(f'Missing chiral center atom {key}')
            named_neighbors=sorted(b if a==name else a for a,b in edges if a==name or b==name)
            coordinates={};substitutions={}
            for neighbor in named_neighbors:
                neighbor_key=(chain,resnum,resname,neighbor)
                if neighbor_key in by_key:
                    coordinates[neighbor]=np.array(by_key[neighbor_key]['xyz_A']);continue
                external=[topo_atoms[index] for index in neighbors[by_topology_key[key]]
                          if (topo_atoms[index]['segment'],topo_atoms[index]['resnum']) != (chain,resnum)
                          and topo_atoms[index]['mass']>1.1]
                if elements[neighbor]=='O' and len(external)==1:
                    external_key=tuple(external[0][k] for k in ['segment','resnum','resname','atom_name'])
                    if external_key in by_key:
                        coordinates[neighbor]=np.array(by_key[external_key]['xyz_A']);substitutions[neighbor]=list(external_key)
            if len(named_neighbors) not in [3,4] or set(named_neighbors)-set(coordinates):
                centers.append(dict(residue=[chain,resnum,resname],atom=name,status='unresolved',neighbors=named_neighbors));continue
            def volume(values,center):
                vectors=[values[n]-center for n in named_neighbors]
                return float(np.linalg.det(np.array(vectors if len(vectors)==3 else [vectors[k]-vectors[3] for k in range(3)])))
            wanted=volume(ideal,ideal[name]);observed=volume(coordinates,np.array(by_key[key]['xyz_A']))
            centers.append(dict(residue=[chain,resnum,resname],atom=name,ccd_config=cc['pdbx_stereo_config'][i],neighbors=named_neighbors,substitutions=substitutions,ideal_signed_volume=wanted,observed_signed_volume=observed,status='consistent' if wanted*observed>0 else 'inverted'))
    return dict(pdb=str(pdb),topology=str(topology),tested_centers=len(centers),
                status='consistent' if centers and all(c['status']=='consistent' for c in centers) else 'failed',
                typed_heavy_atom_count=len(expected_heavy_keys),model_heavy_atom_count=len(observed_heavy_keys),
                centers=centers,actual_interresidue_CNS_bonds=links,topology_total_charge=actual_topology['total_charge'])


def audit_air(path: Path | str) -> dict:
    text=Path(path).read_text()
    assignments=re.findall(r'assign\s*\(resi\s+(\d+)\s+and\s+segid\s+([A-Za-z0-9]+)\)',text,re.I)
    partners=re.findall(r'resi\s+(\d+)\s+and\s+segid\s+([A-Za-z0-9]+)',text,re.I)
    return dict(path=str(path),sha256=sha256(Path(path)),assignment_count=len(assignments),
                assigned_protein_residues=sorted({int(r) for r,s in assignments if s=='A'}),
                glycan_residues=sorted({int(r) for r,s in partners if s=='B'}),
                information_condition='Known protein binding-site AIRs; native-informed, not blind.')


def verify_lock(root: Path | str=ROOT, lock_path: Path | str=DEFAULT_LOCK) -> dict:
    root=Path(root);lock=json.loads((root/lock_path).read_text());checked=[]
    required=['distribution','source_archive','published_attestation','cns_binary','cns_upstream_binary','environment_freeze']
    try:
        records=[lock[key] for key in required]+[lock['license_context'][key] for key in ['haddock_license','cns_disclaimer']]+[lock['official_example'][key] for key in ['upstream_config','execution_config']]
        inputs=lock['official_example']['input_files']
        if len(inputs)!=4:raise HaddockAuditError('Expected four official example input files')
        records+=inputs
        if any(not isinstance(r,dict) or not {'path','sha256','bytes'} <= set(r) for r in records):
            raise HaddockAuditError('Missing mandatory HADDOCK release/input pin fields')
    except (KeyError,TypeError) as exc:
        raise HaddockAuditError('Missing mandatory HADDOCK release/input pin') from exc
    if lock['distribution']['sha256']!=EXPECTED_SDIST_SHA256 or any(lock[key]['sha256']!=EXPECTED_CNS_SHA256 for key in ['cns_binary','cns_upstream_binary']):
        raise HaddockAuditError('Unexpected official release/CNS identity')
    def visit(value):
        if isinstance(value,dict):
            if set(value) & {'path','sha256','bytes'}:
                if not {'path','sha256','bytes'} <= set(value):
                    raise HaddockAuditError('Incomplete source pin record')
                path=root/value['path']
                if not path.is_file() or sha256(path)!=value['sha256'] or path.stat().st_size!=value['bytes']:
                    raise HaddockAuditError(f'HADDOCK release/input checksum mismatch: {value["path"]}')
                checked.append(value['path'])
            for item in value.values():visit(item)
        elif isinstance(value,list):
            for item in value:visit(item)
    visit(lock)
    if lock['version']!='2026.10.1' or lock['source_commit']!='815fcc14100896e35d1795058465cb547f37a8c5':
        raise HaddockAuditError('Unexpected released version/source commit')
    if 'known_binding_site' not in lock['official_example']['information_condition'] or lock['official_example']['sampling']!=20 or lock['execution']['ncores']!=1:
        raise HaddockAuditError('Released example information/sampling/CPU condition changed')
    upstream=(root/lock['official_example']['upstream_config']['path']).read_text()
    executed=(root/lock['official_example']['execution_config']['path']).read_text()
    if executed!=upstream.replace('ncores = 10','ncores = 1'):
        raise HaddockAuditError('Official scientific configuration changed beyond ncores')
    prefix='haddock3-'+lock['source_commit']+'/examples/docking-protein-glycan/'
    with tarfile.open(root/lock['source_archive']['path']) as archive:
        for record in [lock['official_example']['upstream_config']]+inputs:
            name='docking-protein-glycan-test.cfg' if record is lock['official_example']['upstream_config'] else 'data/'+Path(record['path']).name
            payload=archive.extractfile(prefix+name).read()
            if hashlib.sha256(payload).hexdigest()!=record['sha256']:
                raise HaddockAuditError('Example source does not match pinned official archive: '+name)
    return dict(status='verified',checked_paths=len(set(checked)),lock_sha256=sha256(root/lock_path),
                cns_sha256=lock['cns_binary']['sha256'],version=lock['version'])


def workflow_inventory(run_dir: Path | str) -> dict:
    run_dir=Path(run_dir);stages=[]
    for directory in sorted(run_dir.iterdir()):
        if not directory.is_dir() or not re.match(r'^\d+_',directory.name):continue
        pdbs=sorted(set(directory.glob('*.pdb'))|set(directory.glob('*.pdb.gz')))
        stages.append(dict(stage=directory.name,physical_pdb_files=len(pdbs),
            model_inventory=[pdb_inventory(path) for path in pdbs],
            io_present=(directory/'io.json').is_file(),params_present=(directory/'params.cfg').is_file()))
    return dict(run_dir=str(run_dir),stages=stages)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,default=ROOT)
    parser.add_argument('--lock',type=Path,default=DEFAULT_LOCK);parser.add_argument('--run-dir',type=Path)
    args=parser.parse_args(argv)
    try:
        report=verify_lock(args.root,args.lock)
        if args.run_dir:report['workflow']=workflow_inventory(args.run_dir)
    except (ValueError,OSError,KeyError) as exc:
        print(json.dumps(dict(status='failed',error=str(exc))));return 1
    print(json.dumps(report,indent=2));return 0


if __name__=='__main__':raise SystemExit(main())
