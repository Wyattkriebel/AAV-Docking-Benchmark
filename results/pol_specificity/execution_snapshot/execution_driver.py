from pathlib import Path
from datetime import datetime, timezone
import json,os,shutil,sys
root=Path('/Users/alec/AAV-Docking-Benchmark');sys.path.insert(0,str(root))
from aav_benchmark.pol_specificity import prepare_pol_panel,dock_pol_panel,validate_experiment
from aav_benchmark.provenance import sha256,write_json
from aav_benchmark.archive import bundle
config=root/'benchmarks/pol_experiment.json';prepared=root/'data/processed/pol_specificity';output=root/'results/pol_specificity';vina=root/'tools/vina-1.2.7'
expected='802e044255248bdaf5d75675e3105158722c1aed29541bcf38eec01e3494b71f'
def guard():
 assert sha256(config)==expected,'Frozen POL config changed'
 c=json.loads(config.read_text());validate_experiment(c,execution=True);return c
c=guard();receipt=json.loads((output/'freeze_receipt.json').read_text());assert receipt['configuration_sha256']==expected
assert sha256(vina)=='823c2bbacf26d72183861322345f0a89736aca66c8e81054c66f93af5ad623f1'
print('Preparing complete frozen POL panel',flush=True)
panel=prepare_pol_panel(config,prepared)
print('Preparation succeeded; bundling all16 cases',flush=True)
for row in panel['cases']:bundle(prepared/row['prepared_dir'],output/row['case_id'])
shutil.copyfile(prepared/'panel_preparation.json',output/'panel_preparation.json')
shutil.copyfile(config,output/'pol_experiment.json')
files=list(c['implementation_files_sha256'])+['aav_benchmark/archive.py','scripts/run_pol_specificity.py']
snapshot=output/'execution_snapshot';snapshot.mkdir()
record={'schema_version':1,'configuration_sha256':expected,'executable_sha256':sha256(vina),'started_UTC':datetime.now(timezone.utc).isoformat(),'pid':os.getpid(),'workers':1,'cpu_per_run':2,'files':[]}
for name in files:
 source=root/name;destination=snapshot/name;destination.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,destination);record['files'].append({'path':name,'sha256':sha256(source),'snapshot_path':str(destination.relative_to(output))})
shutil.copyfile('/tmp/run_pol_authorized.py',snapshot/'execution_driver.py');record['driver_sha256']=sha256(snapshot/'execution_driver.py')
write_json(output/'execution_snapshot.json',record);guard();print('Starting96 Vina runs, one worker,2CPU each',flush=True)
status={'schema_version':1,'pid':os.getpid(),'configuration_sha256':expected,'status':'running','expected_runs':96,'completed_cases':0,'terminal_runs':0,'successful_runs':0,'failed_runs':0,'started_UTC':record['started_UTC']}
write_json(output/'execution_status.json',status)
for row in panel['cases']:
 guard();case=row['case_id'];dock_pol_panel(config,prepared,output,vina,case_ids=[case])
 runs=json.loads((output/case/'runs.json').read_text());status['completed_cases']+=1;status['terminal_runs']+=len(runs);status['successful_runs']+=sum(r['status']=='ok' for r in runs);status['failed_runs']+=sum(r['status']=='failed' for r in runs);status['last_case']=case;status['updated_UTC']=datetime.now(timezone.utc).isoformat();write_json(output/'execution_status.json',status);print(f"POL progress {status['terminal_runs']}/96: {status['successful_runs']}ok,{status['failed_runs']}failed",flush=True)
guard();status['status']='completed' if status['successful_runs']==96 else 'completed_with_failures';status['finished_UTC']=datetime.now(timezone.utc).isoformat();write_json(output/'execution_status.json',status);print(json.dumps(status),flush=True)
