# AAV Docking Benchmark

Can molecular docking identify experimentally supported glycan sites and poses,
and distinguish molecular recognition from favorable but incorrect interactions?
This repository evaluates those questions separately, with explicit input
information, structural references and assay provenance.

**Status:** executed development benchmarks for pose diagnostics, blind site
localization and assay-backed specificity. The confirmatory study remains open.
The [protocol](PROTOCOL.md) defines the intended study;
[evidence gates](benchmarks/closure_gates.json) record the unfinished scientific
requirements. A successful software run does not establish biological reliability.

## New site and specificity evidence

The [atomic-reference qualification](docs/ATOMIC_REFERENCE_QUALIFICATION.md)
now reviews the seven external development cases against deposited validation
and Fourier maps. Two full modeled ligands and four explicitly mapped cores
pass the operational criteria; one case is rejected. These decisions retain
borderline density and occupancy cautions and do not imply that every scored
atom is uniquely resolved. The frozen old docking results remain unchanged.
Replay the source checks with `python scripts/verify_atomic_references.py`.

The [new structural-source audit](docs/HELDOUT_REFERENCE_QUALIFICATION.md)
retains a fixed twelve-candidate inventory: six full modeled ligands, five
explicit cores and one atomic rejection. A separate
[family and training audit](docs/FAMILY_LEAKAGE_AUDIT.md) identifies substantial
training-family overlap and unknown provenance. These sources support a new
workflow test; they do not establish model-family novelty. The six priority
cases are fixed before docking. Public replay retains all numerical inputs and
recomputes qualification and 56,366 archived sequence comparisons with
`python scripts/verify_public_source_cohort.py`.

The [site-localization report](docs/SITE_LOCALIZATION_RESULTS.md) evaluates
fpocket 4.2.3, P2Rank 2.5.1 and 100 area-weighted surface draws per receptor,
using whole protein inputs without reference-derived crops or search boxes.
Eight published-site complexes provide the primary denominator; independently
solved AAV9 apo is a paired ninth condition. At the frozen 4 Å endpoint, primary
Top-5 recovery is 5/8 for fpocket and 4/8 for P2Rank. Both full-capsid fpocket
executions crashed and remain eligible failures. P2Rank's first matching AAV
site class ranks 11 in holo and 6 in apo. Site references support localization,
with explicit structural cautions; they do not certify atomic poses.

The [Gal9C specificity report](docs/GAL9_SPECIFICITY_RESULTS.md) evaluates six
measured binders against six source-qualified, size/charge-matched C3 epimers
reported as no binding detected at 2 mM. All 36 frozen Vina runs completed.
Mean Top-1 energies give AUROC 0.7222 and favor the measured binder in 5/6
pairs; size and charge each give AUROC 0.5000. This known-site experiment
describes one assay and chemical series. No-binding labels are categorical,
not numerical affinity bounds, and no alternative has atomic pose truth.

The independent [POL specificity experiment](docs/POL_SPECIFICITY_RESULTS.md)
compares natural D-galactose with actually tested D-glucose from one ITC study.
All 96 prespecified runs completed across two anomers, two sites, independently
solved holo/apo receptors and calcium-retained/dry profiles. The primary
equal-anomer mean favors galactose by only 0.0273 kcal/mol; six of seven
sensitivity conditions reverse that preference. The matched anomers also
favor different sugars. This workflow therefore does not robustly distinguish
this measured epimer contrast under the tested preparation choices. Equal
anomer weighting is a model heuristic, and one pair cannot estimate classifier
accuracy or cross-family specificity.

Recompute the new results from archived bytes without installing predictors:

```sh
python scripts/verify_site_predictions.py --output results/site_discovery/prediction_verification.json --check
python scripts/verify_site_experiment.py --numeric-tolerance 1e-10 --surface-area-tolerance 1e-6
python scripts/regenerate_site_report.py --check
python scripts/regenerate_specificity_report.py --check --numeric-tolerance 1e-10
python scripts/regenerate_pol_report.py --check
```

The site verifier restores native sources offline, reparses all 27 predictor
attempts, resamples all 900 surface draws and recomputes geometry and endpoints
in temporary folders. See [execution and replay instructions](docs/NEXT_EXPERIMENT_REPRODUCTION.md),
[reference audit](docs/SITE_REFERENCE_AUDIT.md) and
[independent review](docs/NEXT_PR_INDEPENDENT_REVIEW.md).

## Reproduce the AAV9 diagnostic

Python 3.11 or later, with enough disk space for scientific dependencies and
roughly 150 MB of working inputs. From the repository root:

```sh
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'
python -m pytest -q
python -m aav_benchmark replay aav9 --output results/replay/aav9
```

The replay restores checksum-pinned 7MTZ/7MT0/GAL sources, prepares receptors
and an independently generated ligand, installs the official Vina 1.2.7 binary
for macOS/Linux ARM64 or x86-64, runs three frozen seeds in four conditions,
then writes raw poses, logs, run manifests, pose metrics and summaries under
the requested output folder. Committed evidence lives under
`results/development/aav9/`; use a fresh replay folder to preserve it. Windows
binaries are not currently supported.
Stages can be run separately with `fetch`, `prepare aav9`, `dock aav9`, and
`analyze aav9`. Use `dock aav9 --resume` to reuse only verified completed runs.
Changed inputs or configurations require a new output directory.

This is **true-site local docking**: the ligand supplies the box and receptor
neighborhood. The complete biological assembly is materialized once; the
whole-subunit neighborhood covers the box plus 6 Å. Independently solved apo
and ligand-stripped holo receptors each have complete and monomer conditions.
The monomer is a deliberately incomplete control. Neither seeds nor capsid
copies are independent biological observations.

## Released examples and exploratory assay ranking

```sh
python -m aav_benchmark replay 154L --output results/replay/154l
python -m aav_benchmark replay 1KJL --output results/replay/1KJL
python -m aav_benchmark panel prepare
python -m aav_benchmark panel dock --results-root results/replay
python -m aav_benchmark panel analyze --results-root results/replay
```

The external development inventory is frozen in
[external_cohort.json](benchmarks/external_cohort.json). `replay` also accepts
1GNY, 1QFO, 2J1T, 2VUZ and 1PWB. Unsupported chemistry or reference mismatches
fail visibly; do not remove such cases from the eligible denominator.
Reusing HADDOCK inputs with Vina does not reproduce published HADDOCK results.
See [external preparation](docs/EXTERNAL_BENCHMARK.md) and
[source qualification](docs/EVIDENCE_CURATION.md).

Two additional methods now have executed released-example checks:
[original Vina-Carb 1.0](docs/VINA_CARB_REPRODUCTION.md) and
[HADDOCK3 2026.10.1](docs/HADDOCK_REPRODUCTION.md). Vina-Carb's carbohydrate
terms require preserved ring/residue identities; its ordinary `--score_only`
path omits those terms. A source-linked fixed-coordinate caller verifies active
carbohydrate scoring and within-ligand ranking. HADDOCK completed the official
20-model guided test and retained all preparation changes, including rebuilt
acetyl atoms. These are method validation examples, with supplied site
information, rather than blind structural accuracy results.

```sh
python -m aav_benchmark.additional_method verify-example
python -m aav_benchmark.haddock_replay
```

An independently generated [CCD ideal-ring initialization profile](docs/INPUT_RING_GEOMETRY.md)
now preserves the full graphs, stereochemistry and named ring geometry of seven
glycan inputs while allowing inter-sugar torsions to vary. Earlier input
profiles remain unchanged. The [review amendment](docs/INPUT_RING_REVIEW_AMENDMENT.md)
explicitly records an overwritten derived review and its unrecovered historical
pin; all 32 primary artifacts verify. Run
`python scripts/verify_ring_initialization.py` for read-only saved-product checks.

The Gal3 panel contains three measured binders from one ITC study. It supports
exploratory affinity ranking, with both reducing-end anomers evaluated separately.
It cannot estimate binder/nonbinder classification. The isolated docking domain
differs from the full-length assayed construct, and glycan size already perfectly
ranks this tiny panel. Report that baseline alongside docking scores.

## What the implementation checks

Atom identity and stereochemistry must survive preparation. Every receptor atom
addition, deletion, hydrogen merge and rejected distance-inferred bond is logged.
RMSD preserves pose error in the receptor frame, contacts retain the full holo
denominator, and retained-pose ranking is distinct from internal sampling.
Distance-only clash diagnostics do not certify ring or linkage validity.

Native coordinates, known translations, receptor transformations, missing atoms,
symmetry operators, failures, stale artifacts and missing run denominators have
analytic regression controls. Statistical utilities average seeds within
complexes before families, retain failed eligible cases, and bootstrap paired
whole families. Assay statistics never pool contexts or invent negative labels.

For a portable preparation bundle use `archive aav9`; `restore aav9` restores it
without overwriting changed inputs. Run `analyze` against archived preparations
to regenerate tables from raw poses. `status --require-complete` exits nonzero
while any scientific evidence gate remains unfinished.

Recompute all published structural and assay tables without rerunning docking:

```sh
python scripts/verify_archived_results.py
python scripts/regenerate_report.py --check
```

The verifier restores archived preparation into temporary folders and compares
freshly computed tables byte for byte. The generated
[development report](docs/DEVELOPMENT_RESULTS.md) includes the failure denominator,
native scoring controls, raw-data hashes and limitations. The
[independent replay](docs/INDEPENDENT_REPRODUCTION.md) and
[clean-environment replay](docs/CLEAN_ENVIRONMENT_REPLAY.md) are implementation
checks on the same host. [Method coverage](docs/METHOD_COVERAGE.md) identifies
the comparisons still needed for the full scientific study.

Linux CI explicitly allows `--numeric-tolerance 1e-10` for derived floating-point
rounding across CPU architectures. It still requires exact identities, counts,
decisions and input/pose hashes, and reports the actual numerical differences.
The default development-host verification remains byte-exact.

Site replay has a separate explicit 0.000001 Å² allowance for accumulated
surface-area totals on Linux. This addresses an observed 0.00000004 Å²
cross-platform difference in a roughly 859,673 Å² total. It does not relax
coordinates, distances, candidate identities, scores or success decisions.

Legacy trimmed structures and pocket outputs remain historical artifacts with
unverified coordinate/provenance mappings. They are not blind benchmark evidence.
