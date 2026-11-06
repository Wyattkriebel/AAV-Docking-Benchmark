#!/usr/bin/env python3
"""Reproduce source-only atom correspondence, map sampling, CCD stereo and geometry diagnostics.
No predicted poses, scores or recovery outcomes are read. Run from repository root.
Requires gemmi==0.7.5, numpy==2.4.2; matplotlib only for optional source density slices.
"""
import gzip,hashlib,json,math,itertools,pathlib,xml.etree.ElementTree as ET
import gemmi,numpy as np
OUT=pathlib.Path('data/archive/atomic_reference')
SITE=json.load(open('benchmarks/site_reference_ledger.json'))
POLICY=json.load(open(OUT/'qualification_policy.json'))
def xyz(cat,i):return np.array([float(cat['Cartn_'+a][i]) for a in 'xyz'])
def read_grid(p,kind):
 b=gemmi.as_refln_blocks(gemmi.cif.read_file(str(OUT/f'{p}_validation_{kind}_map_coef.cif.gz')))[0]
 f,phi=('pdbx_FWT','pdbx_PHWT') if kind=='2fo-fc' else ('pdbx_DELFWT','pdbx_DELPHWT')
 g=b.transform_f_phi_to_map(f,phi,sample_rate=3)
 a=np.asarray(g); mean=float(a.mean()); std=float(a.std());a[:]=(a-mean)/std
 return g,dict(shape=list(g.shape),full_cell_mean_before_normalization=mean,full_cell_sigma=std,unit_cell=list(g.unit_cell.parameters),spacegroup=g.spacegroup.hm)
def sample(g,x):return float(g.interpolate_value(gemmi.Position(*map(float,x))))
def ccd(name):
 b=gemmi.cif.read_file(f'data/archive/external/{name}.cif.gz')[0]
 a=b.get_mmcif_category('_chem_comp_atom.'); bo=b.get_mmcif_category('_chem_comp_bond.')
 ideal={n:np.array([float(a['pdbx_model_Cartn_'+c+'_ideal'][i]) for c in 'xyz']) for i,n in enumerate(a['atom_id'])}
 edges=[(bo['atom_id_1'][i],bo['atom_id_2'][i],bo['value_order'][i]) for i in range(len(bo['atom_id_1'])) if not (bo['atom_id_1'][i].startswith('H') or bo['atom_id_2'][i].startswith('H'))]
 return a,ideal,edges
OFFSETS=[np.array(x)*.25 for x in itertools.product(range(-2,3),repeat=3) if sum(v*v for v in x)<=4]
rows=[]
for case in SITE['cases'][1:]:
 p=case['case_id'].lower(); block=gemmi.cif.read_file(str(OUT/f'{p}.cif'))[0]; cat=block.get_mmcif_category('_atom_site.')
 maps={k:read_grid(p,k) for k in ['2fo-fc','fo-fc']}
 xml=ET.fromstring(gzip.open(OUT/f'{p}_validation.xml.gz','rt').read()); xentry=xml.find('Entry').attrib
 conn=block.get_mmcif_category('_struct_conn.'); links=[]
 for i,kind in enumerate(conn.get('conn_type_id',[])):
  if kind!='covale':continue
  left=[conn['ptnr1_auth_asym_id'][i],conn['ptnr1_auth_seq_id'][i],conn['ptnr1_label_atom_id'][i]]
  right=[conn['ptnr2_auth_asym_id'][i],conn['ptnr2_auth_seq_id'][i],conn['ptnr2_label_atom_id'][i]]
  links.append(dict(left=left,right=right,distance_A=float(conn['pdbx_dist_value'][i])))
 # Released source ligand only, never docking input/output.
 released=gemmi.read_pdb_string(gzip.open(f'data/archive/external/{case["case_id"]}_l_b.pdb.gz','rt').read())
 release_atoms={(r.seqid.num,a.name):(r,a) for chain in released[0] for r in chain for a in r if a.element.name not in ['H','D']}
 residues=[]
 for m in case['evaluation_only']['ligand_residue_map']:
  d=m['deposited']; ids=[i for i in range(len(cat['id'])) if cat['auth_asym_id'][i]==d['author_chain'] and cat['auth_seq_id'][i]==d['author_resnum'] and cat['type_symbol'][i] not in ['H','D']]
  coords={cat['label_atom_id'][i]:xyz(cat,i) for i in ids}
  val=next(x for x in xml.findall('ModelledSubgroup') if x.attrib.get('chain')==d['author_chain'] and x.attrib.get('resnum')==d['author_resnum'] and x.attrib.get('resname')==d['resname'])
  aa,ideal,edges=ccd(d['resname']); names=set(coords)
  full={n for i,n in enumerate(aa['atom_id']) if aa['type_symbol'][i] not in ['H','D']}; missing=sorted(full-names)
  stereo=[]
  for i,n in enumerate(aa['atom_id']):
   if aa['pdbx_stereo_config'][i] not in ['R','S'] or n not in names:continue
   ns=sorted([b if a==n else a for a,b,_ in edges if a==n or b==n]); actual={k:coords[k] for k in ns if k in coords}; substitutions={}
   for k in ns:
    if k in actual:continue
    # Glycosidic leaving oxygen is replaced by the external acceptor oxygen.
    candidates=[]
    for l in links:
     for this,other in [(l['left'],l['right']),(l['right'],l['left'])]:
      if this==[d['author_chain'],d['author_resnum'],n]:
       hits=[j for j in range(len(cat['id'])) if [cat['auth_asym_id'][j],cat['auth_seq_id'][j],cat['label_atom_id'][j]]==other]
       if len(hits)==1 and cat['type_symbol'][hits[0]]=='O':candidates.append((other,xyz(cat,hits[0])))
    if len(candidates)==1 and k.startswith('O'):
     actual[k]=candidates[0][1];substitutions[k]=candidates[0][0]
   if len(ns) not in [3,4] or set(ns)-set(actual):
    stereo.append(dict(atom=n,status='unresolved',neighbors=ns));continue
   def det(c,center):
    v=[c[k]-center for k in ns]
    return float(np.linalg.det(np.array(v if len(v)==3 else [v[0]-v[3],v[1]-v[3],v[2]-v[3]])))
   vi=det(ideal,ideal[n]); va=det(actual,coords[n]);stereo.append(dict(atom=n,ccd_config=aa['pdbx_stereo_config'][i],neighbors=ns,substitutions=substitutions,ideal_signed_volume=vi,observed_signed_volume=va,status='consistent' if vi*va>0 else 'inverted'))
  amap=[]
  for i in ids:
   name=cat['label_atom_id'][i];pos=xyz(cat,i);rr,ra=release_atoms[(m['prepared_resnum'],name)]
   density=dict(two_mFo_DFc_atom_sigma=sample(maps['2fo-fc'][0],pos),two_mFo_DFc_max_within_0p5A_sigma=max(sample(maps['2fo-fc'][0],pos+o) for o in OFFSETS),mFo_DFc_atom_sigma=sample(maps['fo-fc'][0],pos))
   amap.append(dict(prepared_chain=m['prepared_chain'],prepared_resnum=m['prepared_resnum'],prepared_resname=m['prepared_resname'],atom_name=name,element=cat['type_symbol'][i],deposited_author_chain=d['author_chain'],deposited_author_resnum=d['author_resnum'],deposited_label_chain=cat['label_asym_id'][i],deposited_label_resnum=cat['label_seq_id'][i],deposited_atom_site_id=cat['id'][i],xyz_A=pos.tolist(),occupancy=float(cat['occupancy'][i]),B_iso_A2=float(cat['B_iso_or_equiv'][i]),released_coordinate_delta_A=float(np.linalg.norm(pos-np.array([ra.pos.x,ra.pos.y,ra.pos.z]))),**density))
  v={k:float(val.attrib[k]) for k in ['rscc','rsr','avgoccu','owab','EDIAm','OPIA'] if k in val.attrib}; flags=[dict(kind=x.tag,**x.attrib) for x in val]
  supported=sum(a['two_mFo_DFc_max_within_0p5A_sigma']>=1 for a in amap); frac=supported/len(amap);neg=min(a['mFo_DFc_atom_sigma'] for a in amap)
  reasons=[]
  for key,op,bound in [('rscc','min',.8),('rsr','max',.25),('avgoccu','min',.9)]:
   if key not in v or (v[key]<bound if op=='min' else v[key]>bound):reasons.append(f'{key}_{op}_gate')
  if frac<.9:reasons.append('2mFo_DFc_coverage_below_90_percent')
  if neg<=-3:reasons.append('negative_difference_density_at_atom')
  if any(x['kind']=='chiral-outlier' for x in flags):reasons.append('wwPDB_chirality_outlier')
  if any(x['kind'] in ['mog-bond-outlier','mog-angle-outlier'] and abs(float(x['Zscore']))>=5 for x in flags):reasons.append('absolute_bond_or_angle_Z_at_least_5')
  if any(x['status']!='consistent' for x in stereo):reasons.append('CCD_signed_volume_stereochemistry_unresolved_or_inverted')
  if any(a['released_coordinate_delta_A']>.000001 for a in amap):reasons.append('released_deposited_coordinate_mismatch')
  residues.append(dict(prepared_resnum=m['prepared_resnum'],prepared_resname=m['prepared_resname'],deposited=d,heavy_atom_count=len(amap),missing_CCD_heavy_atoms=missing,CCD_missing_atom_note='Only named glycosidic leaving oxygen is allowed absent, and its external replacement is checked in stereo/links.',validation=v,validation_flags=flags,stereochemistry=stereo,local_density=dict(supported_atoms=supported,total_atoms=len(amap),fraction_supported=frac,min_atom_difference_sigma=neg),source_gate_pass=not reasons,source_gate_failures=reasons,atom_map=amap,CCD_internal_heavy_bonds=[list(e) for e in edges if e[0] in names and e[1] in names]))
 rows.append(dict(case_id=case['case_id'],maps={k:v[1] for k,v in maps.items()},entry_validation={k:xentry.get(k) for k in ['PDB-resolution','PDB-R','EDS_R','DCC_R','PDB-Rfree','DCC_Rfree','XMLcreationDate','PDB-revision-date']},covalent_links=links,residues=residues))
 print(case['case_id'],[(r['prepared_resnum'],round(r['local_density']['fraction_supported'],3),r['source_gate_failures']) for r in residues])
(OUT/'source_diagnostics.json').write_text(json.dumps(dict(schema_version=1,software=dict(gemmi=gemmi.__version__,numpy=np.__version__),policy=POLICY,local_sphere_sample_count=len(OFFSETS),cases=rows),indent=2)+'\n')
