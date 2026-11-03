"""Source/name-map-only audit. No experimental or actual pose coordinates read."""
from pathlib import Path
import hashlib, json, tarfile, io, tempfile, sys, argparse
import numpy as np
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root', type=Path, required=True, help='Audited repository or public source fixture')
parser.add_argument('--output', type=Path, required=True, help='Fresh JSON path outside the audited root')
args=parser.parse_args()
ROOT=args.root.resolve(strict=True)
OUTPUT=args.output.resolve()
if not ROOT.is_dir():parser.error('--root must be a directory')
if OUTPUT.is_relative_to(ROOT):parser.error('--output must be outside the audited root')
sys.path.insert(0, str(ROOT))
from aav_benchmark.engine import read_poses
from aav_benchmark.paired_method_report import canonical_carb_poses
from aav_benchmark.additional_method import read_carb_poses
from aav_benchmark.paired_common_pose_rescoring import restore_carb_pose_aliases
from aav_benchmark.heldout_pose_report import canonical_mapping
from aav_benchmark.qualified_pose import ParentGeometry
from rdkit import Chem
CASES=('2RDK','3ZWE','4R9F','5HZA','5V6F','8I5U')
def sha(raw):return hashlib.sha256(raw).hexdigest()
def record(path):
 p=ROOT/path;b=p.read_bytes();return {'path':str(path),'sha256':sha(b),'bytes':len(b)}
def chemical_projection(atoms):
 return [{'name':a['name'],'canonical_atom_key':a['canonical_atom_key'],'rdkit_index':a['rdkit_index'],'branch_num':a['branch_num'],'ccd_id':a['ccd_id'],'atom_name':a['atom_name']} for a in atoms]
def normalize(atoms):
 return [{'name':a['name'],'source_key':a.get('source_key',['B',int(a['branch_num']),'',a['ccd_id'],a['atom_name']])} for a in atoms]
source_maps={};source_records=[]
for profile in ('strict','repair'):
 manifest_path=Path(f'results/structural_pose_v3/{profile}_sources_manifest.json')
 manifest=json.loads((ROOT/manifest_path).read_bytes());archive_path=manifest_path.parent/manifest['archive'];archive_raw=(ROOT/archive_path).read_bytes()
 assert sha(archive_raw)==manifest['archive_sha256']
 with tarfile.open(fileobj=io.BytesIO(archive_raw),mode='r:gz') as tar:
  for case in CASES:
   name=f'{case}/ligand_atom_map.json';raw=tar.extractfile(name).read();atoms=json.loads(raw)['atoms']
   expected=next(f for f in manifest['files'] if f['path']==name)
   assert sha(raw)==expected['sha256']
   names=[a['name'] for a in atoms];keys=[a['canonical_atom_key'] for a in atoms];indices=[a['rdkit_index'] for a in atoms]
   assert len(set(names))==len(set(keys))==len(set(indices))==len(atoms)
   assert all(a['canonical_atom_key']==f"{a['branch_num']}:{a['ccd_id']}:{a['atom_name']}" for a in atoms)
   # Read graph fields only: atom symbols/bond indices, never coordinate columns.
   sdf_name=f'{case}/ligand.sdf';sdf_raw=tar.extractfile(sdf_name).read()
   sdf_expected=next(f for f in manifest['files'] if f['path']==sdf_name)
   assert sha(sdf_raw)==sdf_expected['sha256']
   lines=sdf_raw.decode().splitlines();assert 'V2000' in lines[3]
   atom_count=int(lines[3][:3]);bond_count=int(lines[3][3:6])
   elements=[line[31:34].strip() for line in lines[4:4+atom_count]]
   heavy=[i for i,e in enumerate(elements) if e not in ('H','D')]
   assert set(indices)==set(heavy)
   assert all(elements[a['rdkit_index']]==a['element'] for a in atoms)
   bonds=[(int(line[:3])-1,int(line[3:6])-1,int(line[6:9])) for line in lines[4+atom_count:4+atom_count+bond_count]]
   assert all(0<=i<atom_count and 0<=j<atom_count for i,j,_ in bonds)
   source_maps[profile,case]=atoms
   source_records.append({'profile':profile,'case_id':case,'parent_atoms':len(atoms),'source_map_sha256':sha(raw),
    'names':names,'canonical_keys':keys,'rdkit_indices':indices,'names_key_index_bijection':True,
    'metadata_projection_sha256':sha(json.dumps(chemical_projection(atoms),sort_keys=True).encode()),
    'sdf_sha256':sha(sdf_raw),'source_sdf_graph_heavy_index_element_bijection':True,
    'heavy_index_order_from_map':indices,'rdkit_heavy_order_graph_indices':heavy,
    'source_archive':record(archive_path),'source_manifest':record(manifest_path)})
for case in CASES:assert source_maps['strict',case]==source_maps['repair',case]
inputs_path=Path('results/structural_pose_v3/vina_carb_prepared_inputs.json');inputs=json.loads((ROOT/inputs_path).read_bytes());prepared=[];unsupported=[]
for row in inputs['conditions']:
 key=(row['receptor_profile'],row['case_id']);atoms=source_maps[key];names=[a['name'] for a in atoms]
 if row['status']!='prepared':unsupported.append({'condition_id':row['condition_id'],'preparation_status':row['status'],'no_pose_order':True});continue
 p=Path(row['canonical_atom_map']['path']);raw=(ROOT/p).read_bytes();assert sha(raw)==row['canonical_atom_map']['sha256'];actual=json.loads(raw)
 actual=actual.get('atoms',actual.get('canonical_atom_map',actual.get('atom_map'))) if isinstance(actual,dict) else actual
 assert normalize(actual)==normalize(atoms)
 refpath=Path(row['modern_vina_prepared_directory'])/'reference.json';ref=json.loads((ROOT/refpath).read_bytes())
 assert ref['names']==names and normalize(ref['canonical_atom_map'])==normalize(atoms)
 assert ref.get('pose_reference_eligible') is False and 'coordinates' not in ref
 prepared.append({'condition_id':row['condition_id'],'parent_atoms':len(names),'prepared_reference_names_equal_source_order':True,
  'prepared_canonical_map_equal_source_order':True,'canonical_map':record(p),'prepared_reference':record(refpath)})
# Every coordinate below is synthetic and unrelated to any input/native/output.
coords=np.array([[i+0.125, (i*i)%7+0.25, (i*i*i)%11+0.375] for i in range(5)])
# Source-name order deliberately differs from RDKit heavy order and raw record order.
rdkit_order=[3,1,4,0,2];names=[f'L{i:03}' for i in range(5)]
atoms=[{'name':name,'rdkit_index':rdkit_order[i],'canonical_atom_key':f'1:SYN:X{rdkit_order[i]}','branch_num':1,'ccd_id':'SYN','atom_name':f'X{rdkit_order[i]}'} for i,name in enumerate(names)]
synthetic_case={'inventory_entry':{'expected_heavy_atom_count':5},'scored_atom_count':3,
 'scored_atom_map':[{'canonical_atom_key':atoms[i]['canonical_atom_key']} for i in (2,4,0)],
 'qualified_copies':[{'copy_id':'synthetic','scored_atom_count':3,'scored_atom_map':[{'canonical_atom_key':atoms[i]['canonical_atom_key'],'xyz_A':[j+1.,j+2.,j+3.]} for j,i in enumerate((2,4,0))]}]}
reference={'names':names,'canonical_atom_map':atoms};mapping=canonical_mapping(synthetic_case,reference);assert mapping['indices']==[2,4,0]
# Synthetic graph, synthetic conformer; no embedding or modeling.
mol=Chem.MolFromSmiles('CCNCO');conf=Chem.Conformer(5)
for i,p in enumerate(coords):conf.SetAtomPosition(i,p.tolist())
mol.AddConformer(conf);sdf=Chem.MolToMolBlock(mol).encode()+b'\n$$$$\n'
converted={'names':names,'atom_map':[{**a,'source_key':['L',1,'','SYN',a['atom_name']]} for a in atoms]}
graph=ParentGeometry(sdf,converted,mapping)
assert graph.name_to_heavy==rdkit_order
assert np.array_equal(graph.heavy_to_name,np.argsort(rdkit_order))
assert mapping['indices'] in graph.scored_maps
name_coords=coords[rdkit_order]
assert np.allclose(name_coords[graph.heavy_to_name],graph.input_xyz,atol=1e-3,rtol=0)
serials=[17,6,42,8,20]
export={'atoms':[{'serial':serials[i],'generated_name':names[i],'carb_name':f'C{i+1}','carb_resname':'SYN','carb_resnum':1,'hydrogen':False,'atom_type':'C','charge':0.0} for i in range(5)]}
def atomline(serial,name,point,resname='LIG'):
 return f'HETATM{serial:5d} {name:<4s} {resname:>3s} L{1:4d}    {point[0]:8.3f}{point[1]:8.3f}{point[2]:8.3f}  1.00  0.00    {0.0:6.3f} C'
record_order=[4,1,3,0,2]
with tempfile.TemporaryDirectory(prefix='name-order-synthetic-') as temp:
 temp=Path(temp);modern=temp/'modern.pdbqt';carb=temp/'carb.pdbqt';prepared_ligand=temp/'prepared.pdbqt';restored=temp/'restored.pdbqt'
 modern.write_text('\n'.join(['MODEL 1','REMARK VINA RESULT: -1.0 0 0']+[atomline(serials[i],names[i],name_coords[i]) for i in record_order]+['ENDMDL','']))
 carb.write_text('\n'.join(['MODEL 1','REMARK VINA RESULT: -1.0 0 0']+[atomline(serials[i],f'C{i+1}',name_coords[i],'SYN') for i in reversed(record_order)]+['ENDMDL','']))
 prepared_ligand.write_text('\n'.join([atomline(serials[i],names[i],name_coords[i]) for i in range(5)]+['']))
 m=read_poses(modern,names)[0]['coordinates'];assert np.array_equal(m,name_coords)
 decoded=read_carb_poses(carb,export)[0];assert decoded['names']!=names
 c=canonical_carb_poses(carb,export,names)[0]['coordinates'];assert np.array_equal(c,name_coords)
 restore_carb_pose_aliases(prepared_ligand,carb,export,restored)
 restored_parsed=read_poses(restored,names)[0]['coordinates'];assert np.array_equal(restored_parsed,name_coords)
 assert np.array_equal(m,c) and np.array_equal(m[mapping['indices']],name_coords[[2,4,0]])
 assert np.array_equal(c[graph.heavy_to_name],coords)
source_paths=['aav_benchmark/paired_method_report.py','aav_benchmark/heldout_pose_report.py','aav_benchmark/haddock_pose_report.py','aav_benchmark/assembly_control_report.py','aav_benchmark/paired_common_pose_report.py','aav_benchmark/paired_common_pose_rescoring.py','aav_benchmark/qualified_pose.py','aav_benchmark/engine.py','aav_benchmark/metrics.py','aav_benchmark/additional_method.py','aav_benchmark/prospective_vina_carb.py']
result={'schema_version':1,'scope':'source, prepared name/map metadata and synthetic coordinates only','actual_pose_or_native_coordinates_read':False,'metric_or_accuracy_functions_called':False,'native_execution':False,'repository_files_changed':False,
 'source_input_manifest':record(inputs_path),'source_maps':source_records,'prepared_conditions':prepared,'nonprepared_conditions':unsupported,
 'counts':{'source_profile_case_maps':len(source_records),'input_conditions':len(inputs['conditions']),'prepared_name_contracts':len(prepared),'nonprepared_retained':len(unsupported)},
 'synthetic_controls':{'nontrivial_rdkit_to_source_name_order':True,'nontrivial_raw_modern_record_order':True,'modern_explicit_name_reorder':True,'nontrivial_vc_serial_order':True,'vc_explicit_name_reorder':True,'canonical_mask_order':[2,4,0],'parent_geometry_graph_order_roundtrip':True,'identity_automorphism_in_source_name_mask_order':True,'union_vc_alias_serial_roundtrip':True,'no_rmsd_fit_contact_or_geometry_diagnostic_call':True},
 'source_code_pins':[record(p) for p in source_paths]}
out=OUTPUT;out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':'source_name_order_contracts_pass','counts':result['counts'],'synthetic_controls':result['synthetic_controls'],'output':str(out),'sha256':sha(out.read_bytes())},indent=2))
