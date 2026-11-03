from pathlib import Path
import hashlib,json,sys,os
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from aav_benchmark.engine import dock
from aav_benchmark.prospective_vina_carb import validate_frozen_plan,execute_one
from aav_benchmark.provenance import write_json
OUT=ROOT/'results/structural_pose_v3';MODERN=ROOT/'benchmarks/structural_pose_v3_experiment.json';VC=ROOT/'benchmarks/structural_test_vina_carb.json';VF=OUT/'vina_carb_freeze_receipt.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
plan,conditions=validate_frozen_plan(ROOT,VC,VF)
frozen=json.loads((OUT/'freeze_receipt.json').read_text());assert frozen['configuration_sha256']==sha(MODERN)
rows=[]
def progress(status):
 write_json(OUT/'paired_execution_progress.json',{'status':status,'driver_pid':os.getpid(),'planned_eligible_method_seed_slots':216,'planned_native_search_slots_before_CHI_preflight':114,'native_cpu_per_search':2,'concurrent_searches':1,'terminal_eligible_seed_slots':sum(3 if r['method']=='vina_1_2_7' else 1 for r in rows),'atomic_accuracy_evaluated':False,'records':rows})
progress('running')
for condition in conditions:
 for rel,digest in frozen['helper_sources_sha256'].items():
  if sha(ROOT/rel)!=digest:raise ValueError('Frozen helper changed: '+rel)
 assert sha(ROOT/condition['source_condition']['path'])==condition['source_condition']['sha256']
 identity=condition['condition_id'];print('START_MODERN',identity,flush=True)
 path=ROOT/condition['modern_runs'];receipt_path=path.parents[1]/'execution.json'
 if condition['status']=='prepared':
  for key in ['ligand','receptor','preparation']:
   if sha(ROOT/condition[key]['path'])!=condition[key]['sha256']:raise ValueError('Frozen input changed')
  runs=dock(ROOT/condition['modern_vina_prepared_directory'],path.parent,ROOT/'tools/vina-1.2.7',MODERN)
  result={'condition_id':identity,'configuration_sha256':sha(MODERN),'status':'complete','native_search_started':True,'terminal_seed_runs':len(runs),'seed_statuses':{str(r['seed']):r['status'] for r in runs},'source_condition':condition['source_condition'],'runs_path':condition['modern_runs'],'runs_sha256':sha(path)}
 else:
  result={'condition_id':identity,'configuration_sha256':sha(MODERN),'status':'prepared_input_failure','native_search_started':False,'eligible_seed_slots':plan['seeds'],'failure_reason':condition['failure_reason'],'source_condition':condition['source_condition'],'no_runs_or_poses_invented':True}
 result['frozen_input_manifest_sha256']=plan['prepared_inputs']['sha256'];write_json(receipt_path,result)
 rows.append({'condition_id':identity,'method':'vina_1_2_7','status':result['status'],'execution_receipt':str(receipt_path.relative_to(ROOT))});progress('running');print('DONE_MODERN',identity,result['status'],flush=True)
 for seed in plan['seeds']:
  print('START_VC',identity,seed,flush=True)
  result=execute_one(root=ROOT,plan_path=VC,receipt_path=VF,condition_id=identity,seed=seed,installation=ROOT/'tools/vina-carb-1.0-replay-build',output=ROOT/condition['vina_carb_outputs'][str(seed)])
  rows.append({'condition_id':identity,'seed':seed,'method':'vina_carb_1_0','status':result['status'],'execution_receipt':str(Path(condition['vina_carb_outputs'][str(seed)])/'execution.json')});progress('running');print('DONE_VC',identity,seed,result['status'],flush=True)
progress('complete')
