from pathlib import Path
import hashlib,json,sys,os
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from aav_benchmark.engine import dock
from aav_benchmark.provenance import write_json
OUT=ROOT/'results/structural_assembly_pose_v3';CONFIG=ROOT/'benchmarks/structural_assembly_pose_v3_experiment.json';FREEZE=OUT/'freeze_receipt.json';INPUT=OUT/'prepared_conditions.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
frozen=json.loads(FREEZE.read_text());assert frozen['configuration_sha256']==sha(CONFIG)
inputs=json.loads(INPUT.read_text());assert inputs['status']=='complete' and len(inputs['conditions'])==12
rows=[]
def progress(status):write_json(OUT/'execution_progress.json',{'status':status,'driver_pid':os.getpid(),'planned_eligible_seed_slots':36,'planned_native_seed_runs':27,'terminal_eligible_seed_slots':3*len(rows),'actual_native_cpu_per_search':2,'concurrent_control_searches':1,'atomic_or_site_accuracy_evaluated':False,'records':rows})
progress('running')
for row in inputs['conditions']:
 for rel,digest in frozen['helper_sources_sha256'].items():
  if sha(ROOT/rel)!=digest:raise ValueError('Frozen control helper changed')
 condition_path=ROOT/row['condition_path'];assert sha(condition_path)==row['condition_sha256']
 condition=json.loads(condition_path.read_text());identity='/'.join([row['receptor_profile'],row['case_id'],row['method'],'box01']);base=condition_path.parent.parent/'vina';target=base/'box01';print('START',identity,flush=True)
 if row['prepared_jobs']:
  assert len(condition['jobs'])==1 and condition['jobs'][0]['status']=='prepared'
  job=condition['jobs'][0];prepared=condition_path.parent/job['prepared_directory']
  runs=dock(prepared,target,ROOT/'tools/vina-1.2.7',CONFIG)
  receipt={'status':'complete','condition_id':identity,'native_search_started':True,'terminal_seed_runs':len(runs),'seed_statuses':{str(r['seed']):r['status'] for r in runs},'runs_path':str((target/'runs.json').relative_to(ROOT)),'runs_sha256':sha(target/'runs.json')}
 else:
  receipt={'status':'prepared_input_failure','condition_id':identity,'native_search_started':False,'eligible_seed_slots':[104729,130363,155921],'no_runs_or_poses_invented':True,'failure_reason':condition['status']}
 receipt.update(configuration_sha256=sha(CONFIG),source_condition={'path':row['condition_path'],'sha256':row['condition_sha256']},prepared_conditions_sha256=sha(INPUT))
 write_json(base/'execution.json',receipt);rows.append({'condition_id':identity,'status':receipt['status'],'execution_receipt':str((base/'execution.json').relative_to(ROOT))});progress('running');print('DONE',identity,receipt['status'],flush=True)
progress('complete')
