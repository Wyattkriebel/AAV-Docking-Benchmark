from pathlib import Path
import json,hashlib
root=Path('/tmp/aav-qualified-v2-replay-b8a11f1');sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
report_path=root/'results/development_pose_report_release/qualified_pose_v2_summary.json';r=json.loads(report_path.read_text())
dpath=Path('/tmp/aav-old156-independent-numerical-review.json');d=json.loads(dpath.read_text())
npath=Path('/tmp/aav-old-v2-independent-acceptance-numeric.json');n=json.loads(npath.read_text())
assert sha(dpath)=='0ae68d2bfb838ef513ca7670c0620dedbfed979b7df10c498283268d8d20e593'
assert r['denominator']==d['denominator']==n['denominator']
prior={(c['case_id'],c['method']):c for c in d['conditions']};pools=perms=fields=0
for c in n['conditions']:
    for budget,seeds in c['budgets'].items():
        old={s['seed']:s for s in prior[c['case_id'],c['method']]['budgets'][budget]}
        for s in seeds:
            assert s['returned_poses']==old[s['seed']]['returned_poses']
            assert s['canonical_named_diagnostic']==old[s['seed']]['canonical_named_diagnostic'];pools+=1
eligible=[c for c in r['conditions'] if c['atomic_pose_eligible']];seen={};cross=0
for c in eligible:
    for p in c['atom_order_correction']:
        o=next(o for o in d['order_checks'] if (o['case_id'],o['method'],o['box_id'])==(c['complex_id'],c['method'],p['box_id']))
        assert p['raw_box_names']==o['box_name_order']
        assert p['canonical_to_raw_indices']==[p['raw_box_names'].index(x) for x in p['canonical_parent_names']]
        for f in [p['raw_reference'],p['canonical_reference']]:assert sha(root/f['archive'])==f['archive_sha256']
        b=p['geometry_diagnostic_baseline'];assert b==c['geometry_diagnostic_baseline'] and b['exact_executed_box_PDBQT_baseline'] is False
        assert b['source']=='archived_independent_legacy_SDF' and sha(root/b['archive'])==b['archive_sha256'];perms+=1
    for budget,b in c['site_budgets'].items():
        for s in b['seed_evaluations']:
            for p in s['poses']:
                key=c['complex_id'],c['method'],s['seed'],p['box_id'],p['rank'];value=p['global_RMSD_A'],p['local_interface_RMSD_A']
                if key in seen:assert seen[key]==value;cross+=1
                seen[key]=value
for scope in ['full_ligand','modeled_core']:
    for method in ['fpocket','p2rank','oracle']:
        cases=[c for c in eligible if c['atomic_scope']==scope and c['method']==method]
        assert len(cases)==(2 if scope=='full_ligand' else 4)
        for budget in ['1','5']:
            for fit in ['global','local_interface']:
                for rank in ['top1','top5','best_returned']:
                    for threshold in [1,2,3]:
                        f=f'{rank}_success_{threshold}A';values=[sum(s[fit][f] for s in c['site_budgets'][budget]['seed_evaluations'])/3 for c in cases]
                        assert all(c['site_budgets'][budget][fit][f]==v for c,v in zip(cases,values))
                        assert abs(r['aggregates'][scope][method][budget][fit][f]-sum(values)/len(values))<1e-15;fields+=1
assert pools==108 and perms==33 and len(seen)==1920 and cross==800
assert len(r['paired_contrasts'])==48 and all(p['mean_oracle_minus_predicted']==0 and 'no population confidence interval' in p['uncertainty'] for p in r['paired_contrasts'])
for path,p in r['source_sha256'].items():
    f=root/path;assert f.is_file() and f.stat().st_size==p['bytes'] and sha(f)==p['sha256'],path
out={'status':'official_old_v2_independent_acceptance_passed','official_report_sha256':sha(report_path),'prior_diagnostic_sha256':sha(dpath),'numeric_review_sha256':sha(npath),'numeric_driver_sha256':sha('/tmp/aav_old156_v2_acceptance.py'),'contract_driver_sha256':sha(__file__),'raw_models':1920,'pooled_rows':2720,'cross_budget_rows':cross,'per_box_provenance_checks':perms,'prior_diagnostic_seed_pools_exact':pools,'aggregate_fields_independently_checked':fields,'maximum_independent_RMSD_error_A':n['maximum_official_named_RMSD_difference_A'],'source_files_verified':len(r['source_sha256']),'atomic_eligible_pools_per_budget':54,'full_pools_per_budget':18,'core_pools_per_budget':36,'denominator':r['denominator'],'aav_and_154L_null':all(c['atomic_accuracy'] is None for c in r['conditions'] if not c['atomic_pose_eligible']),'1QFO_eligible_zero':all(c['atomic_pose_eligible'] and c['atomic_scope']=='modeled_core' and all(not s['poses'] for b in c['site_budgets'].values() for s in b['seed_evaluations']) for c in r['conditions'] if c['complex_id']=='1QFO'),'actual_report_regenerated_by_reviewer':False,'new_six_case_accuracy_accessed':False,'native_calls':False,'frozen_edits':False}
p=Path('/tmp/aav-old-v2-independent-acceptance.json');p.write_text(json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+'\n');print(json.dumps(out,indent=2))
