#!/usr/bin/env python3
"""Source-only observed-OXT terminal amendment; retains all prior jobs."""
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from aav_benchmark.haddock_cohort_preparation import pin, audit_topology
from aav_benchmark.haddock_method import pdb_atoms
from scripts.preflight_haddock_pairs import prepare_pair, audit_receptor, POLICY as BASE_POLICY, OUTPUT as BASE_OUTPUT

POLICY = ROOT / 'benchmarks/haddock_terminal_amendment.json'
OUTPUT = BASE_OUTPUT / 'observed_termini_v2'


def observed_terminal_rule(atoms):
    chains = {a['chain'] for a in atoms}
    present = {a['chain'] for a in atoms if a['atom_name'] == 'OXT'}
    if not present:
        return False
    if present != chains:
        raise ValueError('Mixed observed OXT chain inventory')
    for chain in chains:
        oxygens = [a for a in atoms if a['chain'] == chain and a['atom_name'] == 'OXT']
        if len(oxygens) != 1 or int(oxygens[0]['resnum']) != max(int(a['resnum']) for a in atoms if a['chain'] == chain):
            raise ValueError('Duplicate or nonterminal observed OXT')
    return True


def main():
    policy = json.loads(POLICY.read_text())
    base_policy = json.loads(BASE_POLICY.read_text())
    jobs, records = [], []
    for profile in policy['profiles_retained']:
        for case in policy['cases_retained']:
            source = ROOT / 'data/processed/heldout_pose_v3' / profile / case / 'receptor.pdb'
            charged = observed_terminal_rule(pdb_atoms(source))
            if not charged:
                original = BASE_OUTPUT / profile / case / 'paired_topology_audit.json'
                records.append(dict(profile=profile, case_id=case, charged_cter=False,
                                    execution='unchanged_source_configuration_reuses_v1_execution',
                                    original_audit=pin(original), audit=json.loads(original.read_text())))
                continue
            destination, prepared = prepare_pair(case, profile, OUTPUT)
            config = destination / 'topology.cfg'
            config.write_text(config.read_text() + '[topoaa.mol1]\ncharged_cter = true\n')
            prepared['config'] = pin(config)
            prepared['terminal_amendment'] = pin(POLICY)
            (destination / 'preparation.json').write_text(json.dumps(prepared, indent=2) + '\n')
            jobs.append((profile, case, destination, prepared))
    OUTPUT.mkdir(parents=True, exist_ok=True)
    plan = dict(schema_version=1, written_before_first_CNS_job=True, source_only=True, no_docking=True,
                policy=pin(POLICY), driver=pin(Path(__file__)),
                jobs=[dict(profile=p, case_id=c, config=prep['config'], preparation=pin(d / 'preparation.json'))
                      for p, c, d, prep in jobs], reused=[r['original_audit'] for r in records])
    (OUTPUT / 'execution_plan.json').write_text(json.dumps(plan, indent=2) + '\n')
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1', MKL_NUM_THREADS='1')
    for profile, case, destination, prepared in jobs:
        start = time.monotonic()
        with (destination / 'stdout.log').open('w') as stdout, (destination / 'stderr.log').open('w') as stderr:
            process = subprocess.run([str(ROOT / 'tools/haddock3-2026.10.1/env/bin/haddock3'), 'topology.cfg'],
                                     cwd=destination, env=environment, stdout=stdout, stderr=stderr, check=False)
        audit = dict(profile=profile, case_id=case, source_only=True, no_docking_performed=True,
                     process_return_code=process.returncode, elapsed_seconds=time.monotonic() - start,
                     glycan=audit_topology(case, destination), receptor=audit_receptor(destination, base_policy))
        (destination / 'paired_topology_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
        records.append(dict(profile=profile, case_id=case, charged_cter=True, execution='source_only_terminal_amendment_topoaa',
                            amended_audit=pin(destination / 'paired_topology_audit.json'), audit=audit))
        print(profile, case, process.returncode, audit['glycan']['status'], audit['receptor']['status'], flush=True)
    (OUTPUT / 'summary.json').write_text(json.dumps(dict(schema_version=1, source_only=True, no_docking_performed=True,
                                                       policy=pin(POLICY), driver=pin(Path(__file__)), cases=records), indent=2) + '\n')


if __name__ == '__main__':
    main()
