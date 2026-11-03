#!/usr/bin/env python3
"""Source-only missing receptor contacts/assembly context for selected native atoms."""
import pathlib,json,gzip,hashlib
import gemmi,numpy as np
OUT=pathlib.Path('data/archive/atomic_reference');site=json.load(open('benchmarks/site_reference_ledger.json'));ledger={c['case_id']:c for c in json.load(open('benchmarks/atomic_reference_ledger.json'))['cases']}
rows=[]
for c in site['cases'][1:]:
 p=c['case_id'];b=gemmi.cif.read_file(str(OUT/f'{p.lower()}.cif'))[0];a=b.get_mmcif_category('_atom_site.');s=gemmi.read_pdb_string(gzip.open(f'data/archive/external/{p}_r_b.pdb.gz','rt').read());rxyz=np.array([[aa.pos.x,aa.pos.y,aa.pos.z] for ch in s[0] for r in ch for aa in r if aa.element.name not in ['H','D']]);lig=[];core=[]
 for i in range(len(a['id'])):
  if a['type_symbol'][i] in ['H','D']:continue
  for m in c['evaluation_only']['ligand_residue_map']:
   if a['auth_asym_id'][i]==m['deposited']['author_chain'] and a['auth_seq_id'][i]==m['deposited']['author_resnum']:
    x=np.array([float(a['Cartn_'+v][i]) for v in 'xyz']);lig.append(x)
    if m['prepared_resnum'] in ledger[p]['scored_prepared_residue_numbers']:core.append(x)
 lig=np.array(lig);core=np.array(core);chain_dist={};missing=[]
 for i in range(len(a['id'])):
  if a['type_symbol'][i] in ['H','D'] or a['label_seq_id'][i] in [False,None,'.','?']:continue
  x=np.array([float(a['Cartn_'+v][i]) for v in 'xyz']);d=float(np.min(np.linalg.norm(lig-x,axis=1)));dc=float(np.min(np.linalg.norm(core-x,axis=1))) if len(core) else None;ch=a['auth_asym_id'][i];chain_dist[ch]=min(chain_dist.get(ch,1000),d)
  if d<4 and float(np.min(np.linalg.norm(rxyz-x,axis=1)))>.001:missing.append(dict(author_chain=ch,author_resnum=a['auth_seq_id'][i],resname=a['label_comp_id'][i],atom_name=a['label_atom_id'][i],alt_id=a['label_alt_id'][i],occupancy=float(a['occupancy'][i]),distance_full_ligand_A=d,distance_scored_core_A=dc))
 row=dict(case_id=p,protein_chain_min_distance_to_full_selected_ligand_A=chain_dist,protein_atoms_absent_from_released_coordinates_within_4A=missing,assembly_generation=b.get_mmcif_category('_pdbx_struct_assembly_gen.'))
 if p=='2VUZ':
  st=gemmi.make_structure_from_block(b);op=list(gemmi.find_spacegroup_by_name(st.spacegroup_hm).operations())[3];xyz=[]
  for ch in st[0]:
   if ch.name!='A':continue
   for r in ch:
    for aa in r:
     f=st.cell.fractionalize(aa.pos);f=op.apply_to_xyz([f.x,f.y,f.z]);x=st.cell.orthogonalize(gemmi.Fractional(*f));xyz.append([x.x,x.y,x.z])
  row['published_dimer_operator_4_555']=dict(operator=op.triplet(),distance_partner_protein_to_scored_core_A=float(np.min(np.linalg.norm(np.array(xyz)[:,None,:]-core[None,:,:],axis=2))))
 rows.append(row)
(OUT/'receptor_context_diagnostics.json').write_text(json.dumps(dict(schema_version=1,source_only=True,contact_cutoff_A=4,cases=rows),indent=2)+'\n')
print('written receptor context diagnostics')
