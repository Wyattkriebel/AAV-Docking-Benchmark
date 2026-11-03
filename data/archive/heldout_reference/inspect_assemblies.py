#!/usr/bin/env python3
"""Source-only assembly and receptor-context inventory; no prediction inputs/outcomes."""
import pathlib,json
import gemmi,numpy as np
OUT=pathlib.Path('data/archive/heldout_reference');j=json.load(open(OUT/'audit_inputs.json'));rows=[]
for c in j['candidates']:
 st=gemmi.read_structure(c['sources']['coordinates']['archive_path']);b=gemmi.cif.read_file(c['sources']['coordinates']['archive_path'])[0];a=b.get_mmcif_category('_atom_site.');selected=c['selected_label_asym_id'];liglabels={selected}
 for e in c['extra_covalent_edges']:
  if e['ptnr1_label_asym_id']==selected:liglabels.add(e['ptnr2_label_asym_id'])
  elif e['ptnr2_label_asym_id']==selected:liglabels.add(e['ptnr1_label_asym_id'])
 ligidx=[i for i in range(len(a['id'])) if a['label_asym_id'][i] in liglabels and a['type_symbol'][i] not in ['H','D']]
 lig=np.array([[float(a['Cartn_'+v][i]) for v in 'xyz'] for i in ligidx]);proteins={};extras={}
 for i in range(len(a['id'])):
  if a['type_symbol'][i] in ['H','D']:continue
  x=np.array([float(a['Cartn_'+v][i]) for v in 'xyz']);d=float(np.min(np.linalg.norm(lig-x,axis=1)));lab=a['label_asym_id'][i]
  if a['label_seq_id'][i] not in [None,False,'.','?']:proteins[lab]=min(proteins.get(lab,1000),d)
  elif lab not in liglabels:
   key=(lab,a['auth_asym_id'][i],a['auth_seq_id'][i],a['label_comp_id'][i]);extras[key]=min(extras.get(key,1000),d)
 assembly_ids=[g['assembly_id'] for g in c['assembly_generators'] if selected in g['asym_id_list'].split(',')];assembly_id=sorted(set(assembly_ids))[0] if assembly_ids else None
 assembly=next((aa for aa in st.assemblies if aa.name==assembly_id),None);generated=[]
 if assembly:
  for g in assembly.generators:
   for op in g.operators:
    generated.append(dict(operator_name=op.name,subchains=list(g.subchains),matrix=op.transform.mat.tolist(),vector=op.transform.vec.tolist()))
 ions=[dict(label_chain=k[0],author_chain=k[1],author_resnum=k[2],resname=k[3],min_distance_selected_full_ligand_A=v,retained=(k[3] in ['CA','MN'] or (c['pdb_id']=='2WNB' and k[3]=='C5P'))) for k,v in extras.items() if k[3]!='HOH']
 row=dict(case_id=c['pdb_id'],selected_assembly_id=assembly_id,selected_ligand_label_chains=sorted(liglabels),assembly_metadata=c['assembly_metadata'],assembly_generators=c['assembly_generators'],expanded_assembly_operators=generated,deposited_protein_label_chain_min_distance_A=proteins,native_nonwater_nonpolymer_components=ions,receptor_policy=dict(protein='All protein subchains and operators of selected deposited biological assembly, including deposited MSE and covalent modifications; no native-neighborhood crop.',retained_cofactors=[x for x in ions if x['retained']],omitted_components='Strip selected glycan/aglycone, other carbohydrate copies, waters, GOL, EDO, sulfate, iodide, TRS and other crystallization additives. Retain all deposited CA/MN; for2WNB retain C5P(CMP) as explicit primary co-ligand state.',alternate_protein_conformers='Generic highest-occupancy conformer, tie lexical altID with blank/common atoms retained; record exact mapping. No ligand-dependent conformer selection.',assembly_scope='Complete deposited selected construct assembly, not full physiological protein/virion/environment.'),deposited_source_organism=dict(entity_src_gen=b.get_mmcif_category('_entity_src_gen.'),entity_src_nat=b.get_mmcif_category('_entity_src_nat.')))
 rows.append(row)
 print(c['pdb_id'],assembly_id,'protein',proteins,'retain',[x for x in ions if x['retained']])
(OUT/'assembly_context.json').write_text(json.dumps(dict(schema_version=1,source_only=True,prediction_outcomes_consulted=False,cases=rows),indent=2)+'\n')
