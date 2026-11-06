"""Separate, source-only chemical-ideal pyranose conformer initialization.

No receptor, native ligand coordinates, docking output or accuracy result is read.
The old independent preparation and every frozen study remain unchanged.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import gzip
import hashlib
import json
import math
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem, DistanceGeometry
from rdkit.Chem import AllChem, rdDistGeom, rdMolTransforms

from .ring_geometry import RING_ATOM_ORDER, compare_to_ideal, cremer_pople6

ROOT=Path(__file__).resolve().parents[1]
POLICY_PATH=Path('benchmarks/ideal_ring_conformer_policy.json')


class RingConformerError(ValueError):
    pass


def _sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def _volume(center,neighbors):
    vectors=np.asarray(neighbors)-center
    return float(np.linalg.det(vectors if len(vectors)==3 else vectors[:3]-vectors[3]))
def _angle_delta(a,b):return abs((a-b+180)%360-180)


def graph_fingerprint(molecule):
    heavy=Chem.RemoveHs(Chem.Mol(molecule));Chem.AssignStereochemistry(heavy,cleanIt=True,force=True)
    atoms=sorted((a.GetProp('canonical_atom_key'),a.GetSymbol(),a.GetFormalCharge(),a.GetIsAromatic()) for a in heavy.GetAtoms())
    bonds=sorted((tuple(sorted((b.GetBeginAtom().GetProp('canonical_atom_key'),b.GetEndAtom().GetProp('canonical_atom_key')))),str(b.GetBondType()),b.GetIsAromatic()) for b in heavy.GetBonds())
    return json.dumps(dict(atoms=atoms,bonds=bonds,isomeric_smiles=Chem.MolToSmiles(heavy,isomericSmiles=True)),sort_keys=True,separators=(',',':'))


def case_schema(case_id,root=ROOT):
    """Filter only chemical fields from frozen metadata, never read source 3D."""
    root=Path(root)
    if case_id=='1GNY':
        path=root/'benchmarks/external_cohort.json';j=json.loads(path.read_text());c=next(c for c in j['cases'] if c['case_id']==case_id)
        residues=[dict(branch_num=str(r['resnum']),ccd_id=r['comp_id']) for r in c['residues']]
        sources={comp:j['ccd_sources'][comp] for comp in sorted({r['ccd_id'] for r in residues})}
        links=c['links']
    else:
        path=root/'benchmarks/heldout_candidate_inventory.json';j=json.loads(path.read_text());c=next(c for c in j['candidates'] if c['pdb_id']==case_id)
        if c['extra_covalent_edges']:raise RingConformerError('External covalent parent requires a separate chemical schema')
        residues=[dict(branch_num=r['branch_num'],ccd_id=r['ccd_id']) for r in c['residues']]
        sources={r['filename'].split('.')[0]:r for r in j['ccd_sources'] if r['filename'].split('.')[0] in {r['ccd_id'] for r in residues}}
        links=c['branch_links']
    return dict(case_id=case_id,residues=residues,branch_links=links,ccd_sources=sources,
                source_metadata=dict(path=str(path.relative_to(root)),sha256=_sha(path)),experimental_coordinates_read=False)


def _component(source,root):
    path=Path(root)/(source.get('archive') or source['archive_path']);raw=path.read_bytes();payload=gzip.decompress(raw) if path.suffix=='.gz' else raw
    if source.get('archive_sha256') and hashlib.sha256(raw).hexdigest()!=source['archive_sha256']:raise RingConformerError('CCD archive checksum mismatch')
    if hashlib.sha256(payload).hexdigest()!=source['sha256']:raise RingConformerError('CCD source checksum mismatch')
    block=gemmi.cif.read_string(payload.decode())[0];a=block.get_mmcif_category('_chem_comp_atom.');b=block.get_mmcif_category('_chem_comp_bond.')
    atoms={name:dict(element=a['type_symbol'][i],charge=int(a['charge'][i]),stereo=a['pdbx_stereo_config'][i],aromatic=a.get('pdbx_aromatic_flag',['N']*len(a['atom_id']))[i]=='Y',ideal=np.array([float(a['pdbx_model_Cartn_'+v+'_ideal'][i]) for v in 'xyz'])) for i,name in enumerate(a['atom_id']) if a['type_symbol'][i] not in ['H','D']}
    bonds=[(b['atom_id_1'][i],b['atom_id_2'][i],b['value_order'][i].upper()) for i in range(len(b['atom_id_1'])) if b['atom_id_1'][i] in atoms and b['atom_id_2'][i] in atoms]
    edit=Chem.RWMol();indices={};orders={'SING':Chem.BondType.SINGLE,'DOUB':Chem.BondType.DOUBLE,'TRIP':Chem.BondType.TRIPLE,'AROM':Chem.BondType.AROMATIC}
    for name,d in atoms.items():
        atom=Chem.Atom(d['element']);atom.SetFormalCharge(d['charge']);atom.SetIsAromatic(d['aromatic']);indices[name]=edit.AddAtom(atom)
    for l,r,o in bonds:edit.AddBond(indices[l],indices[r],orders[o])
    mol=edit.GetMol();Chem.SanitizeMol(mol);conf=Chem.Conformer(mol.GetNumAtoms())
    for name,i in indices.items():conf.SetAtomPosition(i,atoms[name]['ideal'])
    mol.AddConformer(conf);Chem.AssignAtomChiralTagsFromStructure(mol,replaceExistingTags=True);Chem.AssignStereochemistry(mol,cleanIt=True,force=True)
    if not all(name in atoms for name in RING_ATOM_ORDER):raise RingConformerError('Unsupported or incomplete named pyranose ring')
    ring_indices=[indices[n] for n in RING_ATOM_ORDER];ring_sets=[set(r) for r in mol.GetRingInfo().AtomRings()]
    if ring_sets!=[set(ring_indices)]:raise RingConformerError('Only one unambiguous six-membered pyranose per CCD is supported')
    xyz=np.array([atoms[n]['ideal'] for n in RING_ATOM_ORDER]);torsions=[rdMolTransforms.GetDihedralDeg(conf,*[ring_indices[(k+q)%6] for q in range(4)]) for k in range(6)]
    return dict(atoms=atoms,bonds=bonds,mol=mol,indices=indices,ring_xyz=xyz,ring_torsions=torsions,source=dict(path=str(path),archive_sha256=hashlib.sha256(raw).hexdigest(),uncompressed_sha256=hashlib.sha256(payload).hexdigest(),coordinate_fields='pdbx_model_Cartn_{x,y,z}_ideal only'))


def build_parent_graph(schema,root=ROOT):
    """Construct complete mapped chemistry; return zero-conformer parent graph."""
    components={name:_component(source,root) for name,source in schema['ccd_sources'].items()};residues={r['branch_num']:r['ccd_id'] for r in schema['residues']};replacements={};links=[]
    for row in schema['branch_links']:
        if row['value_order'].upper()!='SING':raise RingConformerError('Unsupported branch bond order')
        ends=[f"{row[f'entity_branch_list_num_{i}']}:{row[f'comp_id_{i}']}:{row[f'atom_id_{i}']}" for i in [1,2]];links.append(ends)
        for i in [1,2]:
            key=f"{row[f'entity_branch_list_num_{i}']}:{row[f'comp_id_{i}']}:{row[f'leaving_atom_id_{i}']}"
            if row[f'leaving_atom_id_{i}'] in components[row[f'comp_id_{i}']]['atoms']:replacements[key]=ends[2-i]
    edit=Chem.RWMol();indices={};mapping=[]
    for number,comp in residues.items():
        c=components[comp]
        for name,definition in c['atoms'].items():
            key=f'{number}:{comp}:{name}'
            if key in replacements:continue
            atom=Chem.Atom(c['mol'].GetAtomWithIdx(c['indices'][name]));atom.SetProp('canonical_atom_key',key);atom.SetProp('ccd_atom_name',name)
            i=edit.AddAtom(atom);indices[key]=i;mapping.append(dict(rdkit_index=i,canonical_atom_key=key,branch_num=number,ccd_id=comp,atom_name=name,element=definition['element'],formal_charge=definition['charge'],ccd_stereo_config=definition['stereo']))
        for bond in c['mol'].GetBonds():
            names=list(c['atoms']);l=f'{number}:{comp}:{names[bond.GetBeginAtomIdx()]}';r=f'{number}:{comp}:{names[bond.GetEndAtomIdx()]}'
            if l in indices and r in indices:edit.AddBond(indices[l],indices[r],bond.GetBondType())
    for l,r in links:edit.AddBond(indices[l],indices[r],Chem.BondType.SINGLE)
    mol=edit.GetMol();reverse={i:k for k,i in indices.items()};stereo=[];rings=[]
    for number,comp in residues.items():
        c=components[comp];names=list(c['atoms'])
        for name,definition in c['atoms'].items():
            key=f'{number}:{comp}:{name}'
            if key not in indices or definition['stereo'] not in ['R','S']:continue
            original=c['mol'].GetAtomWithIdx(c['indices'][name]);old_neighbors=[f'{number}:{comp}:{names[a.GetIdx()]}' for a in original.GetNeighbors()];desired=[replacements.get(k,k) for k in old_neighbors];actual=[reverse[a.GetIdx()] for a in mol.GetAtomWithIdx(indices[key]).GetNeighbors()]
            if set(actual)!=set(desired) or len(actual) not in [3,4]:raise RingConformerError('CCD stereo neighborhood changed unexpectedly')
            permutation=[desired.index(k) for k in actual];inversions=sum(permutation[i]>permutation[j] for i in range(len(permutation)) for j in range(i+1,len(permutation)))
            if inversions%2:mol.GetAtomWithIdx(indices[key]).InvertChirality()
            neighbors=sorted(old_neighbors);volume=_volume(definition['ideal'],[c['atoms'][k.split(':')[2]]['ideal'] for k in neighbors])
            stereo.append(dict(center=key,neighbors=[replacements.get(k,k) for k in neighbors],ideal_signed_volume=volume,ccd_stereo_config=definition['stereo']))
        rings.append(dict(branch_num=number,ccd_id=comp,indices=[indices[f'{number}:{comp}:{name}'] for name in RING_ATOM_ORDER],ideal_xyz=c['ring_xyz'].tolist(),ideal_torsions=c['ring_torsions'],ideal_pucker=cremer_pople6(c['ring_xyz'])))
    Chem.SanitizeMol(mol);Chem.AssignStereochemistry(mol,cleanIt=True,force=True)
    if len(Chem.GetMolFrags(mol))!=1 or any(s=='?' for _,s in Chem.FindMolChiralCenters(mol,includeUnassigned=True,useLegacyImplementation=False)):raise RingConformerError('Disconnected parent or unassigned stereochemistry')
    mol.RemoveAllConformers();mol=Chem.AddHs(mol);mol.SetProp('expected_chemical_graph',graph_fingerprint(mol))
    provenance=dict(isomeric_smiles=Chem.MolToSmiles(Chem.RemoveHs(mol),isomericSmiles=True),heavy_atom_count=len(mapping),formal_charge=Chem.GetFormalCharge(mol),declared_links=links,leaving_replacements=replacements,ccd_sources={k:v['source'] for k,v in components.items()},experimental_coordinates_read=False,input_graph_conformer_count=mol.GetNumConformers())
    return mol,mapping,rings,stereo,provenance


def validate_conformer(mol,conf_id,mapping,rings,stereo,policy):
    gates=policy['validity_gates'];xyz=mol.GetConformer(conf_id).GetPositions();keys={a['canonical_atom_key']:a['rdkit_index'] for a in mapping};failures=[];stereo_rows=[];ring_rows=[]
    if not mol.HasProp('expected_chemical_graph') or graph_fingerprint(mol)!=mol.GetProp('expected_chemical_graph'):failures.append('full_parent_bond_charge_or_isomeric_graph_mismatch')
    if len(mapping)!=mol.GetNumHeavyAtoms() or set(keys)!=set(a.GetProp('canonical_atom_key') for a in mol.GetAtoms() if a.GetAtomicNum()!=1):failures.append('full_parent_atom_map_mismatch')
    for s in stereo:
        v=_volume(xyz[keys[s['center']]],[xyz[keys[k]] for k in s['neighbors']]);passed=math.isfinite(v) and v*s['ideal_signed_volume']>0 and abs(v)>=gates['min_abs_named_stereocenter_volume_A3'] and abs(v)>=gates['min_named_stereocenter_volume_fraction_of_CCD']*abs(s['ideal_signed_volume']);stereo_rows.append(dict(**s,observed_signed_volume=v,status='consistent' if passed else 'inverted_or_degenerate'))
        if not passed:failures.append('full_parent_stereochemistry')
    for ring in rings:
        points=xyz[ring['indices']];ideal=np.array(ring['ideal_xyz']);cp=cremer_pople6(points);comparison=compare_to_ideal(points,ideal);pair=float(np.max(np.abs(np.linalg.norm(points[:,None]-points[None,:],axis=2)-np.linalg.norm(ideal[:,None]-ideal[None,:],axis=2))));torsions=[rdMolTransforms.GetDihedralDeg(mol.GetConformer(conf_id),*[ring['indices'][(k+q)%6] for q in range(4)]) for k in range(6)];td=max(_angle_delta(a,b) for a,b in zip(torsions,ring['ideal_torsions']));bad=[]
        checks=[(comparison['pucker_vector_angle_to_CCD_ideal_deg'],'max_full_named_pucker_vector_angle_to_CCD_deg'),(comparison['theta_difference_to_CCD_ideal_deg'],'max_named_theta_difference_to_CCD_deg'),(cp['chair_pole_distance_deg'],'max_nearest_chair_pole_distance_deg'),(abs(comparison['Q_difference_to_CCD_ideal_A']),'max_Q_difference_to_CCD_A'),(comparison['max_abs_ring_bond_difference_to_CCD_ideal_A'],'max_ring_bond_difference_to_CCD_A'),(comparison['max_abs_ring_angle_difference_to_CCD_ideal_deg'],'max_ring_angle_difference_to_CCD_deg'),(pair,'max_ring_pair_distance_difference_to_CCD_A'),(td,'max_ring_torsion_difference_to_CCD_deg')]
        for value,label in checks:
            if value is None or not math.isfinite(value) or value>gates[label]:bad.append(label)
        for value,label in [(cp['ring_bond_lengths_A'],'ring_bond_absolute_range_A'),(cp['ring_angles_deg'],'ring_angle_absolute_range_deg')]:
            if min(value)<gates[label][0] or max(value)>gates[label][1]:bad.append(label)
        if bad:failures.extend('ring'+ring['branch_num']+':'+b for b in bad)
        ring_rows.append(dict(branch_num=ring['branch_num'],ccd_id=ring['ccd_id'],pucker=cp,comparison=comparison,max_pair_distance_difference_A=pair,max_torsion_difference_deg=td,failures=bad))
    distance=Chem.GetDistanceMatrix(mol);table=Chem.GetPeriodicTable();worst=1000.0
    for i in range(len(mapping)):
        for j in range(i):
            if distance[i,j]<=gates['nonbonded_graph_distance_exclusion']:continue
            ratio=float(np.linalg.norm(xyz[i]-xyz[j]))/(table.GetRvdw(mol.GetAtomWithIdx(i).GetAtomicNum())+table.GetRvdw(mol.GetAtomWithIdx(j).GetAtomicNum()));worst=min(worst,ratio)
    if worst<gates['minimum_nonbonded_heavy_separation_vdw_sum_fraction']:failures.append('nonbonded_heavy_clash')
    return dict(valid=not failures,failures=sorted(set(failures)),rings=ring_rows,stereochemistry=stereo_rows,min_nonbonded_heavy_vdw_sum_ratio=worst)


def generate_conformer(molecule,mapping,rings,stereo,policy):
    """Bounded independent ETKDG/MMFF search; never uses supplied coordinates."""
    molecule=Chem.Mol(molecule);molecule.RemoveAllConformers();sample=policy['sampling'];constrain=policy['constraints']
    if not AllChem.MMFFHasAllMoleculeParams(molecule):raise RingConformerError('Full parent unsupported by MMFF94')
    props=AllChem.MMFFGetMoleculeProperties(molecule,mmffVariant=constrain['mmff_variant']);attempts=[];best=None
    for slot in range(sample['candidate_slots']):
        candidate=Chem.Mol(molecule);params=AllChem.ETKDGv3();params.randomSeed=sample['seed_base']+slot;params.numThreads=1;params.maxIterations=sample['max_embedding_iterations_per_slot'];params.enforceChirality=True;params.useSmallRingTorsions=sample['use_small_ring_torsions'];params.useRandomCoords=sample['use_random_coordinates'];params.trackFailures=True
        code=AllChem.EmbedMolecule(candidate,params);record=dict(slot=slot,seed=params.randomSeed,embedding_return_code=code,embedding_failure_counts=list(params.GetFailureCounts()))
        if code<0:record['status']='embedding_failed';attempts.append(record);continue
        ff=AllChem.MMFFGetMoleculeForceField(candidate,props,confId=0)
        for ring in rings:
            ideal=np.array(ring['ideal_xyz'])
            for i in range(6):
                for j in range(i):
                    d=float(np.linalg.norm(ideal[i]-ideal[j]));w=constrain['ring_pair_bound_halfwidth_A'];ff.MMFFAddDistanceConstraint(ring['indices'][i],ring['indices'][j],False,d-w,d+w,constrain['distance_force_constant'])
            for k,target in enumerate(ring['ideal_torsions']):ff.MMFFAddTorsionConstraint(*[ring['indices'][(k+q)%6] for q in range(4)],False,target-constrain['ring_torsion_halfwidth_deg'],target+constrain['ring_torsion_halfwidth_deg'],constrain['torsion_force_constant'])
        ff.Initialize();status=ff.Minimize(maxIts=constrain['max_mmff_iterations'],forceTol=constrain['force_tolerance'],energyTol=constrain['energy_tolerance']);validation=validate_conformer(candidate,0,mapping,rings,stereo,policy);record.update(mmff_return_code=status,constrained_mmff_energy=ff.CalcEnergy(),validation=validation)
        if status!=0:validation['failures'].append('mmff_not_converged');validation['valid']=False
        energy=float(AllChem.MMFFGetMoleculeForceField(candidate,props,confId=0).CalcEnergy());record['unconstrained_mmff_energy_at_selected_geometry']=energy;record['status']='valid' if validation['valid'] and math.isfinite(energy) else 'invalid';attempts.append(record)
        if record['status']=='valid' and (best is None or (energy,slot)<best[:2]):best=(energy,slot,Chem.Mol(candidate))
    report=dict(candidate_slots=sample['candidate_slots'],attempts=attempts,valid_candidate_count=sum(a['status']=='valid' for a in attempts),status='selected' if best else 'no_valid_conformer',selected_slot=best[1] if best else None,selected_unconstrained_mmff_energy=best[0] if best else None,coordinates_source='Independent ETKDG with chemical-ideal intraring bounds/dihedrals; no source/native pose or coordMap')
    return (best[2] if best else None),report


def run_profile(root=ROOT,output=None,cases=None):
    root=Path(root);policy=json.loads((root/POLICY_PATH).read_text());output=Path(output or root/'results/ideal_ring_conformer_feasibility'/policy['profile_id']);output.mkdir(parents=True,exist_ok=True);records=[]
    for case_id in cases or policy['cases']:
        destination=output/case_id;destination.mkdir(exist_ok=True)
        if (destination/'report.json').exists():raise RingConformerError('Refuse to overwrite an earlier conformer trial')
        schema=case_schema(case_id,root);mol,mapping,rings,stereo,provenance=build_parent_graph(schema,root);selected,search=generate_conformer(mol,mapping,rings,stereo,policy)
        report=dict(schema_version=1,case_id=case_id,profile_id=policy['profile_id'],policy_sha256=_sha(root/POLICY_PATH),algorithm_source_sha256=_sha(Path(__file__)),ring_math_source_sha256=_sha(root/'aav_benchmark/ring_geometry.py'),schema=schema,chemical_graph=provenance,atom_map=mapping,ring_templates=rings,stereo_tests=stereo,search=search,no_docking_or_accuracy_evaluation=True)
        if selected:
            selected.SetProp('profile_id',policy['profile_id']);selected.SetProp('canonical_atom_map',json.dumps(mapping));Chem.CreateAtomStringPropertyList(selected,'canonical_atom_key');writer=Chem.SDWriter(str(destination/'ligand.sdf'));writer.write(selected);writer.close();reread=Chem.SDMolSupplier(str(destination/'ligand.sdf'),removeHs=False)[0];revalidation=validate_conformer(reread,0,mapping,rings,stereo,policy);report['serialized_sdf_validation']=revalidation
            if not revalidation['valid']:raise RingConformerError('Serialized selected SDF fails chemical/physical gates')
            report['selected_sdf_sha256']=_sha(destination/'ligand.sdf')
        (destination/'report.json').write_text(json.dumps(report,indent=2)+'\n');records.append(dict(case_id=case_id,status=search['status'],valid_candidate_count=search['valid_candidate_count'],selected_slot=search['selected_slot'],report_sha256=_sha(destination/'report.json')));print(records[-1],flush=True)
    return dict(schema_version=1,profile_id=policy['profile_id'],cases=records,no_docking_or_accuracy_evaluation=True)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=ROOT);p.add_argument('--output',type=Path);p.add_argument('--cases',nargs='+');a=p.parse_args(argv)
    report=run_profile(a.root,a.output,a.cases);print(json.dumps(report,indent=2));return 0 if all(r['status']=='selected' for r in report['cases']) else 1


if __name__=='__main__':raise SystemExit(main())
