# AAV Docking Benchmark

Can molecular docking identify experimentally supported glycan sites and poses,
and distinguish molecular recognition from favorable but incorrect interactions?
This repository evaluates those questions separately, with explicit input
information, structural references and assay provenance.

**Status:** reproducible development diagnostics, not a completed confirmatory
benchmark. The [protocol](PROTOCOL.md) defines the intended study;
[evidence gates](benchmarks/closure_gates.json) record the unfinished scientific
requirements. A successful software run does not establish biological reliability.

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

Legacy trimmed structures and pocket outputs remain historical artifacts with
unverified coordinate/provenance mappings. They are not blind benchmark evidence.
