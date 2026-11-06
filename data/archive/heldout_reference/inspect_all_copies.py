#!/usr/bin/env python3
"""Replay all deposited full chemical copies before outcomes; retain other known sites."""
from pathlib import Path
import json
from aav_benchmark.heldout_reference import audit_all,audit_candidate,derive_copy_entries
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'data/archive/heldout_reference'
j=json.load(open(OUT/'audit_inputs.json'));policy=json.load(open(ROOT/j['policy_path']));rows=[];copyinputs=[]
primary=audit_all(root=ROOT);(OUT/'source_diagnostics.json').write_text(json.dumps(primary,indent=2)+'\n')
for e in j['candidates']:
 derived=derive_copy_entries(e,ROOT);copies=[]
 for ce in derived['same_full_chemical_entries']:
  diag=audit_candidate(ce,ROOT,policy);copies.append(dict(copy_id=ce['copy_id'],assembly_memberships=ce['assembly_memberships'],diagnostics=diag));copyinputs.append(ce)
  print(ce['copy_id'],[(r['branch_num'],r['source_gate_failures']) for r in diag['residues']])
 rows.append(dict(case_id=e['pdb_id'],copies=copies,other_carbohydrate_sites=derived['other_carbohydrate_sites'],assembly_generators=derived['assembly_generators'],assembly_operators=derived['assembly_operators']))
(OUT/'copy_inputs.json').write_text(json.dumps(dict(schema_version=1,source_only=True,prediction_outcomes_consulted=False,copies=copyinputs),indent=2)+'\n')
(OUT/'all_copy_diagnostics.json').write_text(json.dumps(dict(schema_version=1,source_only=True,prediction_outcomes_consulted=False,cases=rows),indent=2)+'\n')
