import json,pathlib,hashlib,gzip,datetime
OUT=pathlib.Path('data/archive/atomic_reference')
source=json.load(open(OUT/'source_diagnostics.json')); site={c['case_id']:c for c in json.load(open('benchmarks/site_reference_ledger.json'))['cases'][1:]};retr=json.load(open(OUT/'retrieval_manifest.json'));manual=json.load(open(OUT/'manual_adjudications.json'))
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
ledger=dict(schema_version=1,curation_date='2026-11-03',status='frozen_source_only_atomic_development_reference_qualification_before_new_pose_metric_review',cohort_denominator=7,cohort_selection='All seven already released external development cases; no replacement or case selection from outcomes.',source_only=True,prediction_outcomes_consulted=False,policy_archive=str(OUT/'qualification_policy.json'),policy_sha256=sha(OUT/'qualification_policy.json'),manual_adjudication_archive=str(OUT/'manual_adjudications.json'),manual_adjudication_sha256=sha(OUT/'manual_adjudications.json'),scope='Conservative operational atomic-reference qualification of selected modeled atoms in existing development systems. These are model-phased experimental references; no independent confirmatory, unseen-family, biological specificity, binding-energy or complete-assembly certification.',source_record='benchmarks/sources/atomic_reference_sources.json',required_policy=source['policy'],score_map_rule='Use exactly scored_atom_map. Same full released parent ligand graph is preserved; core means a partial coordinate endpoint, not a truncated or chemically repaired ligand. Excluded atoms must never be silently scored as atomic truth.',cases=[])
for c in source['cases']:
 p=c['case_id'];old=site[p];adj=manual['decisions'][p];rr=c['residues'];passed={r['prepared_resnum'] for r in rr if r['source_gate_pass'] and str(r['prepared_resnum']) not in adj['additional_residue_blocks']}
 maps={(r['deposited']['author_chain'],r['deposited']['author_resnum']):r['prepared_resnum'] for r in rr}; edges=[]
 for l in c['covalent_links']:
  a=maps.get(tuple(l['left'][:2]));b=maps.get(tuple(l['right'][:2]));
  if a is not None and b is not None:edges.append((a,b))
 comps=[];pending=set(passed)
 while pending:
  comp={min(pending)};queue=list(comp)
  while queue:
   a=queue.pop()
   for l,r in edges:
    b=r if l==a else l if r==a else None
    if b in passed and b not in comp:comp.add(b);queue.append(b)
  pending-=comp;comps.append(sorted(comp))
 counts={r['prepared_resnum']:r['heavy_atom_count'] for r in rr};comps.sort(key=lambda cc:(-sum(counts[r] for r in cc),tuple(cc)))
 full=len(passed)==len(rr) and adj['full_ligand_complete'];core=comps[0] if comps else [];eligible=bool(core);scope='full_ligand' if full else 'modeled_core' if eligible else 'rejected'
 atom_map=[a for r in rr if r['prepared_resnum'] in core for a in r['atom_map']]
 for a in atom_map:a['reference_atom_key']=f"{a['prepared_chain']}:{a['prepared_resnum']}:{a['prepared_resname']}:{a['atom_name']}"
 archives={r['kind']:dict(url=r['url'],archive=r['archive'],sha256=r['sha256'],bytes=r['bytes']) for r in retr if r['pdb_id']==p and r['status']=='retrieved'}
 ligpath=f'data/archive/external/{p}_l_b.pdb.gz';recpath=f'data/archive/external/{p}_r_b.pdb.gz'
 case=dict(case_id=p,pdb_id=p,role='development_atomic_reference' if eligible else 'development_source_rejected_atomic_reference',atomic_pose_eligible=eligible,full_ligand_eligible=full,modeled_core_eligible=eligible and not full,qualification=scope,full_ligand_atomic_reference_status='accepted_with_cautions' if full else 'rejected_or_inconclusive',site_only_reference_remains_available=True,source_family_label=old['source_family_label'],qualified_claim='Atomic pose of exactly the declared modeled atom subset' if eligible else 'No experimentally qualified atomic pose in this audit',decision_reason=adj['reason'],scored_prepared_residue_numbers=core,scored_atom_count=len(atom_map),released_ligand_heavy_atom_count=sum(counts.values()),excluded_prepared_residue_numbers=sorted(set(counts)-set(core)),parent_ligand_graph_preserved=True,chemical_repair_authorized=False,sources=archives,released_ligand_source=dict(archive=ligpath,compressed_sha256=sha(ligpath),uncompressed_sha256=hashlib.sha256(gzip.open(ligpath,'rb').read()).hexdigest(),upstream_commit='344e4a41a300b68b7ee3585895a460c863545ddb'),released_receptor_source=dict(archive=recpath,compressed_sha256=sha(recpath),uncompressed_sha256=hashlib.sha256(gzip.open(recpath,'rb').read()).hexdigest()),receptor_scope=old['receptor_scope'],map_evidence=c['maps'],entry_validation=c['entry_validation'],covalent_links=c['covalent_links'],residue_qualification=[{k:v for k,v in r.items() if k not in ['atom_map','CCD_internal_heavy_bonds']}|dict(additional_manual_block=adj['additional_residue_blocks'].get(str(r['prepared_resnum'])),scored=r['prepared_resnum'] in core) for r in rr],full_atom_diagnostic_archive=str(OUT/'source_diagnostics.json'),receptor_context_diagnostic_archive=str(OUT/'receptor_context_diagnostics.json'),scored_atom_map=atom_map,local_density_review='Actual wwPDB validation coefficient maps reconstructed and sampled in3D; ring-plane contour figures also rendered and visually reviewed. No ligand omit map was reconstructed.',no_claims=['independent_confirmatory_validation','biological_specificity','complete_biological_assembly','unseen_fold_generalization','unique_atomic_conformation'])
 ledger['cases'].append(case)
ledger['qualification_counts']={k:sum(c['qualification']==k for c in ledger['cases']) for k in ['full_ligand','modeled_core','rejected']}
pathlib.Path('benchmarks/atomic_reference_ledger.json').write_text(json.dumps(ledger,indent=2)+'\n')
print(json.dumps([(c['case_id'],c['qualification'],c['scored_prepared_residue_numbers'],c['scored_atom_count']) for c in ledger['cases']]))
