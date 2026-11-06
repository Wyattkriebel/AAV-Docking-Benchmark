"""Official OLD156 v2 acceptance; same pinned independent named arithmetic, no evaluator imports.

For v2 the expected published order IS canonical; compare named values directly.
This review never generates an official report or invokes a native executable.
"""
from pathlib import Path
import gzip
import hashlib
import json
import math
import re

import gemmi
import numpy as np

ROOT = Path('/tmp/aav-qualified-v2-replay-b8a11f1')
REPORT = ROOT/'results/development_pose_report_release/qualified_pose_v2_summary.json'
OUT = Path('/tmp/aav-old-v2-independent-acceptance-numeric.json')
AA = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
SEEDS = [104729, 130363, 155921]
SOURCES = {}


def read(path):
    path = Path(path)
    raw = path.read_bytes()
    SOURCES[str(path.relative_to(ROOT))] = hashlib.sha256(raw).hexdigest()
    return raw


def obj(path):
    return json.loads(read(path))


def normalized(key):
    return str(key[0]), int(key[1]), str(key[2]).strip(), str(key[3]), str(key[4])


def pdb_atoms(raw, protein=False):
    atoms = {}
    for line in raw.decode().splitlines():
        if line[:6].strip() not in {'ATOM', 'HETATM'}:
            continue
        name, residue = line[12:16].strip(), line[17:20].strip()
        element = line[76:78].strip() or name.lstrip('0123456789')[0]
        if element.upper() in {'H', 'D'} or (protein and residue not in AA):
            continue
        key = (line[21].strip(), int(line[22:26]), line[26].strip(), residue, name)
        assert key not in atoms, key
        atoms[key] = np.array([float(line[a:a+8]) for a in (30,38,46)])
    return atoms


def cif_protein(raw):
    c = gemmi.cif.read_string(raw.decode()).sole_block().get_mmcif_category('_atom_site.')
    out = {}
    for i, element in enumerate(c['type_symbol']):
        assert element.upper() not in {'H','D'}
        key = (c['auth_asym_id'][i],int(c['auth_seq_id'][i]),str(c['pdbx_PDB_ins_code'][i] or '').strip(),c['label_comp_id'][i],c['label_atom_id'][i])
        assert key not in out
        out[key] = np.array([float(c[x][i]) for x in ('Cartn_x','Cartn_y','Cartn_z')])
    return out


def rigid_fit(a,b):
    a,b=np.array(a),np.array(b); ac,bc=a.mean(0),b.mean(0)
    u,_,vt=np.linalg.svd((a-ac).T@(b-bc))
    correction=np.eye(3); correction[2,2]=np.linalg.det(u@vt)
    rotation=u@correction@vt; translation=bc-ac@rotation
    return rotation,translation


def rms(a,b):
    return float(np.sqrt(np.sum((np.asarray(a)-np.asarray(b))**2)/len(a)))


def archived(directory,filename,revision=None):
    manifest=obj(directory/'manifest.json')
    matches=[f for f in manifest['files'] if f['filename']==filename and (revision is None or f.get('revision')==revision)]
    assert len(matches)==1,(directory,filename,revision)
    record=matches[0]; raw=read(directory/record['archive'])
    assert hashlib.sha256(raw).hexdigest()==record['archive_sha256']
    plain=gzip.decompress(raw)
    assert hashlib.sha256(plain).hexdigest()==record['sha256']
    return plain


def models(raw,names):
    answer=[]; current=None
    for line in raw.decode().splitlines():
        if line.startswith('MODEL '):
            assert current is None
            current={'rank':int(line.split()[1]),'atoms':{}}
        elif line.startswith('REMARK VINA RESULT:'):
            current['score']=float(line.split()[3])
        elif line[:6].strip() in {'ATOM','HETATM'}:
            name=line[12:16].strip()
            if line[77:].strip() in {'H','HD','HS'}:
                continue
            assert name in names, 'unexpected heavy identity'
            if name in names:
                assert name not in current['atoms'], 'duplicate heavy identity'
                current['atoms'][name]=[float(line[a:a+8]) for a in (30,38,46)]
        elif line.startswith('ENDMDL'):
            assert set(current['atoms'])==set(names), 'missing heavy identity'
            current['xyz']=np.asarray([current['atoms'][n] for n in names]);answer.append(current);current=None
    assert current is None
    assert [a['rank'] for a in answer]==list(range(1,len(answer)+1))
    return answer


def indicators(v):
    return {'minimum_RMSD_A':min(v) if v else None,**{f'{rank}_{threshold}A':any(x<=threshold for x in v[:n])
        for rank,n in [('top1',1),('top5',5),('best',len(v))] for threshold in (1,2,3)}}


report=obj(REPORT)
assert SOURCES[str(REPORT.relative_to(ROOT))]=='831e889196cd5755943ffef5bca6844e322c89cafffb161cc96159aceadb97d2'
ledger=obj(ROOT/'benchmarks/atomic_reference_ledger.json')
records={c['case_id']:c for c in ledger['cases']}
eligible=[c for c in report['conditions'] if c['atomic_pose_eligible']]
assert len(report['conditions'])==27 and len(eligible)==18
assert sum(c['atomic_scope']=='full_ligand' for c in eligible)==6
assert sum(c['atomic_scope']=='modeled_core' for c in eligible)==12
excluded=[c for c in report['conditions'] if not c['atomic_pose_eligible']]
assert len(excluded)==9 and all(c['atomic_accuracy'] is None for c in excluded)
assert {c['complex_id'] for c in excluded}=={'aav9','154L'}
assert not ROOT.joinpath('tools').exists() and not ROOT.joinpath('data/processed').exists()
contexts={}; frame_checks=[]
for cid in sorted({c['complex_id'] for c in eligible}):
    record=records[cid]; native=np.asarray([a['xyz_A'] for a in record['scored_atom_map']])
    receptor=record['released_receptor_source']; packed=read(ROOT/receptor['archive'])
    assert hashlib.sha256(packed).hexdigest()==receptor['compressed_sha256']
    protein=pdb_atoms(gzip.decompress(packed),protein=True)
    source=cif_protein(archived(ROOT/f'results/site_discovery/{cid}/inputs','receptor.cif','original'))
    assert set(source)==set(protein)
    protein_keys=list(protein); xyz=np.asarray([protein[k] for k in protein_keys])
    contact_mask=np.linalg.norm(native[:,None,:]-xyz[None,:,:],axis=2)<4.
    residues={protein_keys[j][:4] for j in np.where(contact_mask)[1]}
    fits={}; row=next(c for c in eligible if c['complex_id']==cid)
    for label in ('global','local_interface'):
        keys=[k for k in protein if k[-1] in {'N','CA','C','O'} and (label=='global' or k[:4] in residues)]
        a=np.asarray([source[k] for k in keys]); b=np.asarray([protein[k] for k in keys]);rotation,translation=rigid_fit(a,b)
        stored=row['receptor_alignments'][label]
        assert len(keys)==stored['matched_backbone_atom_count']
        assert np.max(np.abs(rotation-np.array(stored['rotation'])))<1e-12
        assert np.max(np.abs(translation-np.array(stored['translation'])))<1e-11
        assert rms(a@rotation+translation,b)<1e-10
        fits[label]=(rotation,translation)
    ctrl=row['native_controls']; assert ctrl['native_contact_count']==int(contact_mask.sum())
    assert rms(native,native)==ctrl['native_RMSD_A']==0
    assert rms(native+[3,4,0],native)==ctrl['known_translation_RMSD_A']==5
    ref=None; names=None; indices=None
    if cid!='1QFO':
        ref=json.loads(archived(ROOT/f'results/development/{cid}/inputs','reference.json'))
        names=ref['names']; mapping={normalized(a['source_key']):a['name'] for a in ref['atom_map']}
        keys=[(a['prepared_chain'],int(a['prepared_resnum']),'',a['prepared_resname'],a['atom_name']) for a in record['scored_atom_map']]
        indices=[names.index(mapping[k]) for k in keys]
        assert indices==row['atom_mapping']['indices']
        assert rms(np.array(ref['coordinates'])[indices],native)<1e-12
        native_parent=pdb_atoms(gzip.decompress(read(ROOT/record['released_ligand_source']['archive'])))
        assert all(np.linalg.norm(np.array(ref['coordinates'])[i]-native_parent[normalized(a['source_key'])])<1e-12 for a in ref['atom_map'] for i in [names.index(a['name'])])
        assert row['chemical_scored_atom_map_count']==1
    contexts[cid]=(native,fits,ref,names,indices)
    frame_checks.append({'case_id':cid,'protein_heavy_atoms':len(protein),'source_native_max_atom_delta_A':float(max(np.linalg.norm(source[k]-protein[k]) for k in source)),
        'fit_backbone_atoms':{k:row['receptor_alignments'][k]['matched_backbone_atom_count'] for k in fits},'native_contact_count':int(contact_mask.sum()),
        'native_control_RMSD_A':0.,'translated_native_control_RMSD_A':5.,'canonical_mask_matches_released_source':cid!='1QFO','full_graph_map_count':row['chemical_scored_atom_map_count']})

conditions=[]; max_error=0.; max_official_reproduction_error=0.; mismatched_pose_rows=0; raw_models=0; pooled_rows=0; adversarial={}; order_checks=[]
for condition in eligible:
    cid,method=condition['complex_id'],condition['method']; native,fits,ref,names,indices=contexts[cid]
    preparation=obj(ROOT/f'results/blind_pose/{cid}/{method}/condition_preparation.json')
    raw={seed:[] for seed in SEEDS}
    for job in preparation['jobs']:
        if job['status']!='prepared':continue
        directory=ROOT/f'results/blind_pose/{cid}/{method}'/job['id']
        blind_reference=json.loads(archived(directory/'inputs','reference.json'))
        blind_names=blind_reference['names']; assert set(blind_names)==set(names)
        native_in_box_order=np.asarray(ref['coordinates'])[[names.index(n) for n in blind_names]]
        order_checks.append({'case_id':cid,'method':method,'box_id':job['id'],'canonical_and_box_atom_order_equal':blind_names==names,'box_name_order':blind_names,
            'native_pose_wrong_order_RMSD_A':rms(native_in_box_order[indices],native),'native_pose_canonical_order_RMSD_A':rms(np.asarray(ref['coordinates'])[indices],native)})
        runs=obj(directory/'runs.json')
        if method=='oracle': assert np.linalg.norm(np.mean(ref['coordinates'],axis=0)-job['center_A'])<1e-12
        for run in runs:
            if run['status']!='ok':continue
            data=read(directory/run['pose_file']); assert hashlib.sha256(data).hexdigest()==run['pose_sha256']
            parsed=models(data,names); raw_models+=len(parsed)
            assert len(parsed)==run['returned_poses']
            for pose in parsed:
                pose['reported_order_xyz']=np.asarray([pose['atoms'][n] for n in blind_names])
                pose.update(box_id=job['id'],candidate_class_rank=job['candidate_class_rank'],candidate_original_rank=job['candidate_original_rank'])
                raw[run['seed']].append(pose)
    result={'case_id':cid,'method':method,'scope':condition['atomic_scope'],'budgets':{}}
    for budget in (1,5):
        seeds=[]
        for seed in SEEDS:
            pool=sorted((p for p in raw[seed] if p['candidate_class_rank']<=budget),key=lambda p:(p['score'],p['candidate_class_rank'],p['rank']))
            declared=next(s for s in condition['site_budgets'][str(budget)]['seed_evaluations'] if s['seed']==seed)
            assert len(pool)==len(declared['poses'])
            values={label:[] for label in fits}; wrong_values={label:[] for label in fits}
            for ordinal,(pose,published) in enumerate(zip(pool,declared['poses']),1):
                assert (pose['box_id'],pose['rank'],pose['candidate_class_rank'])==(published['box_id'],published['rank'],published['candidate_class_rank'])
                assert published['pooled_rank']==ordinal and pose['score']==published['vina_score']
                for label,(rotation,translation) in fits.items():
                    value=rms((pose['xyz']@rotation+translation)[indices],native); values[label].append(value)
                    wrong=value;wrong_values[label].append(wrong)
                    error=abs(value-published[label+'_RMSD_A']);max_error=max(max_error,error)
                    mismatch=abs(wrong-published[label+'_RMSD_A']);max_official_reproduction_error=max(max_official_reproduction_error,mismatch)
                    assert mismatch<1e-11,(cid,method,seed,budget,ordinal,label,wrong,published[label+'_RMSD_A'])
                    mismatched_pose_rows+=int(label=='global' and error>1e-8)
                    assert published[label+'_chemical_atom_map_index']==0
                pooled_rows+=1
            for label in fits:
                v=wrong_values[label]; a=declared[label]
                for rank,n in [('top1',1),('top5',5),('best_returned',len(v))]:
                    for threshold in (1,2,3):assert a[f'{rank}_success_{threshold}A']==any(x<=threshold for x in v[:n])
                assert a['returned_pose_count']==len(v)
                assert a['best_returned_RMSD_A']==(min(v) if v else None) or abs(a['best_returned_RMSD_A']-min(v))<1e-11
            seeds.append({'seed':seed,'returned_poses':len(pool),'published_order':{label:indicators(v) for label,v in wrong_values.items()},'canonical_named_diagnostic':{label:indicators(v) for label,v in values.items()}})
        result['budgets'][str(budget)]=seeds
    conditions.append(result)
    if cid=='1KJL' and method=='oracle':
        p=raw[SEEDS[0]][0]; xyz=p['xyz']; rotation,translation=fits['global']; baseline=rms((xyz@rotation+translation)[indices],native)
        angle=.71; Q=np.array([[math.cos(angle),-math.sin(angle),0],[math.sin(angle),math.cos(angle),0],[0,0,1]]);shift=np.array([21.,-13.,8.])
        # Use the actual independent source protein backbone, never ligand atoms, for this frame-control fit.
        rec=records[cid]['released_receptor_source']; prot=pdb_atoms(gzip.decompress(read(ROOT/rec['archive'])),protein=True)
        backbone=np.asarray([v for k,v in prot.items() if k[-1] in {'N','CA','C','O'}])
        rr,tt=rigid_fit(backbone@Q+shift,backbone)
        invariant=rms(((xyz@Q+shift)@rr+tt)[indices],native)
        assert abs(invariant-baseline)<1e-11
        sample=read(ROOT/f'results/blind_pose/{cid}/{method}/box01/neighborhood_{SEEDS[0]}.pdbqt')
        first=(sample.decode().split('ENDMDL',1)[0]+'ENDMDL\n').splitlines()
        atom_rows=[i for i,line in enumerate(first) if line[:6].strip() in {'ATOM','HETATM'}]
        reversed_lines=list(first)
        for i,line in zip(atom_rows,reversed([first[i] for i in atom_rows])):reversed_lines[i]=line
        reparsed=models(('\n'.join(reversed_lines)+'\n').encode(),names)[0]['xyz']
        assert np.array_equal(reparsed,xyz)
        heavy=next(i for i in atom_rows if first[i][12:16].strip() in names)
        extra=first[heavy][:12]+'Z999'+first[heavy][16:]
        variants={'atom_loss':first[:heavy]+first[heavy+1:], 'duplicate_name':first[:heavy]+[first[heavy]]+first[heavy:],
                  'extra_heavy_name':first[:heavy]+[extra]+first[heavy:]}
        rejected=[]
        for name,lines in variants.items():
            try:models(('\n'.join(lines)+'\n').encode(),names)
            except AssertionError:rejected.append(name)
        assert set(rejected)==set(variants)
        adversarial={'common_frame_rigid_transform_RMSD_delta_A':abs(invariant-baseline),'native_known_translation_detected_A':5.,
            'atom_record_order_reversal_preserves_named_coordinates':True,'rejected_mutations':rejected,'pool_order_checked_before_RMSD':True}

assert len(conditions)==18
assert all(len(c['budgets'][str(b)])==3 for c in conditions for b in (1,5))
assert all(not s['published_order'][label]['best_2A'] for c in conditions for b in (1,5) for s in c['budgets'][str(b)] for label in ('global','local_interface'))
out={'schema_version':1,'status':'official_old_v2_raw_named_RMSDs_pass','released_summary_sha256':SOURCES[str(REPORT.relative_to(ROOT))],
 'released_raw_commit':'b8a11f1af13e6390063b0312456493c4fa94c175','frozen_reader_modified':False,'benchmark_metric_functions_imported':False,
 'new_six_case_accuracy_accessed':False,'native_executables_used':False,'denominator':report['denominator'],'atomic_eligible_seed_pools_per_budget':54,
 'full_seed_pools_per_budget':18,'core_seed_pools_per_budget':36,'raw_atomic_eligible_models_checked':raw_models,'pooled_pose_rows_checked':pooled_rows,
 'maximum_canonical_named_vs_report_RMSD_difference_A':max_error,'maximum_official_named_RMSD_difference_A':max_official_reproduction_error,
 'pooled_rows_with_canonical_atom_correspondence_mismatch':mismatched_pose_rows,'order_checks':order_checks,
 'diagnostic_only_not_a_replacement_frozen_analysis':True,'frame_and_native_checks':frame_checks,'conditions':conditions,'adversarial_controls':adversarial,
 'source_sha256':SOURCES,'review_driver_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
OUT.write_text(json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+'\n')
print(json.dumps({k:out[k] for k in ('status','raw_atomic_eligible_models_checked','pooled_pose_rows_checked','maximum_canonical_named_vs_report_RMSD_difference_A','maximum_official_named_RMSD_difference_A','pooled_rows_with_canonical_atom_correspondence_mismatch','adversarial_controls')},indent=2))
