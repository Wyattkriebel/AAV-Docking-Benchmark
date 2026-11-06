#!/usr/bin/env python3
"""Archive source provenance; rendered QA/replay reports are derived outputs."""
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parents[3];O=ROOT/'data/archive/heldout_reference'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def record(p):return dict(archive=str(p.relative_to(ROOT)),sha256=sha(p),bytes=p.stat().st_size)
old=json.load(open(ROOT/'data/archive/heldout_candidates/source_manifest.json'));files=[]
for r in old['sources']:
 p=ROOT/r['archive_path'];assert sha(p)==r['archive_sha256'];files.append(record(p))
for p in sorted(O.rglob('*')):
 if p.is_file() and 'density_slices' not in p.parts and p not in [O/'source_manifest.json',O/'replay_verification.json']:files.append(record(p))
for path in ['benchmarks/heldout_candidate_inventory.json','data/archive/heldout_candidates/source_manifest.json','data/archive/atomic_reference/qualification_policy.json']:
 p=ROOT/path;files.append(record(p))
unique={r['archive']:r for r in files};ledger=ROOT/'benchmarks/heldout_reference_ledger.json';manifest=dict(schema_version=1,curation_date='2026-11-03',source_only=True,prediction_outcomes_consulted=False,ledger_archive=str(ledger.relative_to(ROOT)),ledger_sha256=sha(ledger),original_candidate_source_record='data/archive/heldout_candidates/source_manifest.json',original_candidate_pin_count=len(old['sources']),new_experimental_retrievals=json.load(open(O/'retrieval_manifest.json')),primary_paper_retrievals=json.load(open(O/'papers/retrieval_manifest.json')),archived_files=list(unique.values()),derived_density_figures=[record(p) for p in sorted((O/'density_slices').glob('*.png'))],license_context=dict(PDB_archive='CC0 1.0 Universal per https://www.wwpdb.org/about/usage-policies; raw coordinates, structure factors, CCDs and validation products preserve original author citations.',original_candidate_release_metadata='Upstream Apache-2.0 and ProtGlycanDock CC-BY-4.0 contexts are preserved in original inventory; source qualification is independent of released benchmark membership.',primary_articles='Separately licensed publisher/PMC archives; full text/abstract/supplement access scopes remain explicit in primary_paper_evidence.json. No blanket CC0 assertion for articles.'),reconstruction='Actual wwPDB weighted coefficient Fourier maps normalized by full-cell sigma, sampling at frozen33point0.5A spheres; no independent re-refinement or newly generated omit maps.',derived_output_rule='Density PNGs and later replay report are outputs rather than experimental replay inputs. Their checksums are recorded separately; report never hashes itself.')
if (O/'papers/primary_paper_evidence.json').exists():manifest['primary_paper_evidence_archive']='data/archive/heldout_reference/papers/primary_paper_evidence.json'
(O/'source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print('source input pins',len(unique),'originalcandidatepins',len(old['sources']))
