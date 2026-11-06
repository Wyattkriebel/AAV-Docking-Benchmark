#!/usr/bin/env python3
"""Read prepared topology/coordinates and AIRs in CNS, without moving atoms."""
import json
import os
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from aav_benchmark.haddock_cohort_preparation import pin, _checked


def main():
    base = ROOT / 'results/haddock_paired_cohort'
    plan_path = base / 'prepared_plan.json'
    plan = json.loads(plan_path.read_text())
    destination = base / 'AIR_syntax_preflight'
    destination.mkdir(exist_ok=True)
    jobs = []
    for condition in plan['conditions']:
        if condition['status'] != 'ready':
            continue
        _checked(ROOT / condition['AIR']['path'], condition['AIR']['sha256'])
        metadata = json.loads((ROOT / condition['AIR_metadata']['path']).read_text())
        topology = (ROOT / condition['source_preflight']['audit']['path']).parent / 'run-topology/0_topoaa'
        stem = condition['condition_id'].replace('/', '_')
        script = 'set echo=false end\n'
        for molecule in ['receptor', 'glycan']:
            script += f'structure @@{topology / (molecule + "_haddock.psf")} end\n'
        for molecule in ['receptor', 'glycan']:
            script += f'coordinates @@{topology / (molecule + "_haddock.pdb")}\n'
        script += (f'noe reset nrestraints=12000000 class ambi @@{ROOT / condition["AIR"]["path"]} end\n'
                   'noe print threshold=999999 end\ndisplay AIR_PARSE_COMPLETED\nstop\n')
        inp, out = destination / (stem + '.inp'), destination / (stem + '.out')
        inp.write_text(script)
        jobs.append((condition, metadata, inp, out))
    (destination / 'execution_plan.json').write_text(json.dumps(dict(schema_version=1, source_only=True,
        written_before_CNS_syntax_jobs=True, topology_or_coordinates_modified=False, no_pose_generation=True,
        source_plan=pin(plan_path), driver=pin(Path(__file__)),
        jobs=[dict(condition_id=c['condition_id'], CNS_input=pin(i), AIR=c['AIR']) for c, m, i, o in jobs]), indent=2) + '\n')
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    records = []
    binary = ROOT / plan['CNS_binary']['path']
    _checked(binary, plan['CNS_binary']['sha256'])
    for condition, metadata, inp, out in jobs:
        with inp.open('rb') as source, out.open('wb') as target:
            process = subprocess.run([str(binary)], stdin=source, stdout=target, stderr=subprocess.STDOUT,
                                     env=environment, check=False)
        text = out.read_text()
        actual = [int(n) for n in re.findall(r'SELRPN:\s+(\d+) atoms have been selected', text)]
        expected = [len(metadata['protein_heavy_atom_keys']), len(metadata['full_glycan_heavy_atom_keys'])]
        passed = process.returncode == 0 and actual == expected and 'AIR_PARSE_COMPLETED' in text and '-ERR' not in text and re.search(r'of\s+1 NOEs', text) is not None
        records.append(dict(condition_id=condition['condition_id'], status='passed' if passed else 'failed',
                            process_return_code=process.returncode, actual_selection_counts=actual,
                            expected_selection_counts=expected, CNS_input=pin(inp), CNS_output=pin(out)))
    (destination / 'summary.json').write_text(json.dumps(dict(schema_version=1, source_only=True,
        no_pose_generation=True, topology_or_coordinates_modified=False, CNS_binary=plan['CNS_binary'],
        conditions=records), indent=2) + '\n')
    print(json.dumps(dict(checked=len(records), passed=sum(r['status'] == 'passed' for r in records)), indent=2))
    return 0 if all(r['status'] == 'passed' for r in records) else 1


if __name__ == '__main__':
    raise SystemExit(main())
