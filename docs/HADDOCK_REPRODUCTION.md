# HADDOCK3 released carbohydrate workflow reproduction

**HADDOCK3 2026.10.1 successfully ran its official protein–glycan test on this
macOS ARM64 machine.** The complete workflow finished with exit code 0 in
4 minutes 45 seconds using one CPU: 20 rigid-body models, 12 models selected for
flexible refinement, and 12 final models in six clusters. This establishes an
operational glycan-capable method, beyond successful installation. It is a
known-interface guided reproduction, not a blind or held-out benchmark result.

The pre-execution lock is
[`benchmarks/haddock.lock.json`](../benchmarks/haddock.lock.json), SHA256
`5eec702d354597aa56b41a2fda30479bd8b481d35ee9498bb2ec0201f8f62bcb`.
The machine-readable observed result is
[`reproduction_summary.json`](../results/haddock_reproduction/reproduction_summary.json).
No frozen reference qualification files were changed, and no new cohort case
was run with HADDOCK.

## Release, installation and CNS provenance

The [official release](https://github.com/haddocking/haddock3/releases/tag/2026.10.1)
resolves to commit `815fcc14100896e35d1795058465cb547f37a8c5`. The downloaded
[PyPI source distribution](https://pypi.org/project/haddock3/2026.10.1/) has SHA256
`2edf1f7efef5d9d76217f0c378c14781461a5f729c56ad4b91b1926f0d15421c`.
The PyPI-published attestation, signing certificate and GitHub tag response were
archived. The attestation subject digest and certificate's publishing workflow /
source commit agree with the expected release; this audit does not claim a
separate local Sigstore certificate-chain/transparency-log verification.

Installation used the isolated environment `tools/haddock3-2026.10.1/env`.
The project `.venv` was not modified. The exact installed dependency versions,
installation log and source archives are retained. The package was built from
the pinned PyPI sdist; the examples came from the pinned GitHub source archive.

The release includes `src/haddock/cns/bin/arm64-darwin.bin`, a native ARM64 CNS
1.3U executable, 4,171,768 bytes, SHA256
`e6febf74eccb740f8bc5d49955cc8fc202025d5d642e129e7738e7280ad95a88`.
The installed executable has the same hash. A `stop` preflight returned 0 in
0.32 seconds before the docking workflow. No CNS source compilation or
registration was needed. The binary uses the existing Homebrew GCC 16.2.0
`libquadmath.0.dylib`; its resolved path/hash and the helper executable hashes
are recorded in [`runtime_audit.json`](../results/haddock_reproduction/runtime_audit.json).

HADDOCK3 is distributed under Apache 2.0. CNS has a distinct usage condition:
the [official installation documentation](https://www.bonvinlab.org/haddock3-user-manual/install-cns.html)
and [pinned disclaimer](https://github.com/haddocking/haddock3/blob/815fcc14100896e35d1795058465cb547f37a8c5/DISCLAIMER.md)
say it is free for nonprofit applications and requires a proper license for
commercial use. No registration, new legal affirmation or claim of commercial
entitlement was made. The older source-compilation instructions requiring a CNS
download form/password were not used.

## Exact official workflow and information supplied

The executed source is the released
[`docking-protein-glycan-test.cfg`](https://github.com/haddocking/haddock3/blob/815fcc14100896e35d1795058465cb547f37a8c5/examples/docking-protein-glycan/docking-protein-glycan-test.cfg).
Its input is unbound 1LMQ receptor plus a GLYCAM-generated tetrasaccharide.
The copied scientific configuration, sampling, restraints, force-field defaults,
clustering, ranking and analysis stages were unchanged. The sole main-workflow
adaptation was `ncores = 10` → `ncores = 1`, coordinated with the other jobs.
BLAS/OpenMP thread environment variables were also capped at one.

| Stage | Actual models / behavior |
|---|---|
| `topoaa` | Prepared protein and glycan, with typed CNS topology |
| `rigidbody` | 20 models, as specified by the official test |
| RMSD clustering and selection | Four clusters; 12 models available/selected across them |
| `flexref` | All 12 selected models refined with the released defaults |
| Final clustering and selection | Six clusters, 12 final models |
| `caprieval` | Released reference-based analyses, including final `allatoms=true` |

All fourteen stages produced their expected I/O and parameter records. The test
uses 20 initial models, not the 1000-model full production workflow; its result
must not be described as a production-quality benchmark comparison.

The [official example documentation](https://github.com/haddocking/haddock3/blob/815fcc14100896e35d1795058465cb547f37a8c5/docs/pages/examples.md)
explicitly describes known binding-site guidance. `ambig.tbl` contains **15 AIR
assignments** from protein residues
35, 46, 52, 57, 58, 59, 62, 63, 75, 101, 102, 103, 107, 108 and 109 to glycan
residues 1–4. The RMSD-clustering residue lists contain the same interface
information. The native target is supplied to `caprieval` for evaluation.
These inputs place the run in an **oracle / known-interface information condition**.
The published score ranking was retained; models were not selected by their
native-fit metric.

## Actual atom preparation and chemical checks

The official input glycan has NDG1 plus NAG2–4: 109 total atoms and **56 heavy
atoms**. Its acetyl atoms use GLYCAM names `C2N`, `CME`, `O2N`, and reducing NDG1
lacks O1. The target uses standard `C7`, `C8`, `O7` and includes O1, totaling
57 glycan heavy atoms.

The receptor retains all 999 heavy-atom names and coordinates exactly.
Actual glycan CNS preparation produces **111 typed atoms / 57 heavy atoms** with
114 physical bonds and a neutral total charge within floating-point roundoff.
The graph is connected, has four sugar rings, and contains the three expected
NAG C1 → preceding sugar O4 glycosidic bonds. The CNS MTF file uses the `.psf`
extension but contains mmCIF; explicit terminal `-1` sentinel records are excluded
from physical atom/bond/charge counts.

The preparation is **not merely an atom rename**:

- All 44 shared standard heavy-atom names retain exactly their input coordinates.
- Twelve GLYCAM acetyl heavy-atom names disappear. Twelve standard acetyl atoms
  and the missing reducing O1 appear.
- Semantic alias-to-standard coordinate displacements range up to **4.539 Å**;
  NAG3 O2N → O7 is the largest. This is measured rebuilding, not preservation of
  the supplied acetyl coordinates.

A separate topology-only diagnostic enabled `debug=true`, `clean=false` to retain
CNS inputs/stdout, without changing scientific topology parameters or executing
docking. It finished in 3.88 seconds. Its 111 prepared atom records are identical
to the main workflow's preparation. The retained CNS output contains 37 coordinate
read errors, including all **12 heavy acetyl aliases**, confirming that they were
not recognized and their standard counterparts were built. See
[`coordinate_read_audit.json`](../results/haddock_reproduction/topology_diagnostic/coordinate_read_audit.json).
This diagnostic is separate from the unchanged official docking reproduction.

Every prepared glycan heavy atom is matched to the actual CNS topology. Named CCD
handedness checks use actual glycosidic bonds for missing leaving-oxygen
substitutions. All **20 chiral centers** are consistent in the prepared glycan
and in all 20 rigid-body, 12 refined and 12 final model records. These are chemical
inventory/handedness checks, not a density qualification of the example's native
reference. Tests also detect nonchiral atom loss, rather than relying only on
chiral-center counts.

Static parameter review finds default declarations for the priority-case sugars
MAN, BMA, BGC, GAL, GLA and FUC, as well as NAG, NDG and A2G; MSE, CA and MN have
declarations in their corresponding topology files. FRU, OPO and C5P do not have
default declarations in this release's scanned topology set. Static presence does
not prove successful preparation of a new complex. No new cohort case or
unreviewed ligand/cofactor parameterization was attempted.

## Sampling and ranking diagnostic

These are the upstream final `13_caprieval` self-evaluation values, retained to
show what this actual guided workflow did:

| Selection | HADDOCK score | Score rank | Ligand RMSD | DockQ |
|---|---:|---:|---:|---:|
| First by published score | −117.275 | 1 | 11.700 Å | 0.215 |
| Best sampled native fit by DockQ, reported retrospectively | −78.264 | 11 | 2.764 Å | 0.722 |

The stronger native-fit sample was not ranked first. This separates generation
of a better pose from ranking it successfully, even in this supplied-interface
example. It does not establish relative performance versus Vina/Vina-Carb,
validated atomic accuracy, biological specificity or generalization. HADDOCK
scores are method-specific weighted scores and must not be compared numerically
to another method's scores as binding energies. A comparable new-cohort restraint
and preparation protocol must be frozen before HADDOCK enters that comparison.

## Reproduction and retained artifacts

The actual main command ran from
`results/haddock_reproduction/official_glycan_test`:

```sh
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1 \
  /Users/alec/AAV-Docking-Benchmark/tools/haddock3-2026.10.1/env/bin/haddock3 \
  reproduce.cfg
```

For another execution, copy the retained `data` and configs into a fresh working
directory so the existing `run1-test` artifacts remain intact. Release archives,
input/config pins and the exact dependency freeze are in the lock/source review.
All effective parameters, PDB/PSF outputs, stage I/O, native evaluation tables and
workflow logs remain under
[`results/haddock_reproduction/`](../results/haddock_reproduction/).
The installed environment and downloaded executables are under `tools`.

Read-only verification and report regeneration:

```sh
.venv/bin/python -m aav_benchmark.haddock_replay
.venv/bin/python -m pytest tests/test_haddock_method.py -q
```

The portable replay verifies all 423 retained artifacts and nine retained lock
pins, reconstructs preparation, topology, model handedness, stage inventories
and the final raw evaluation table, and compares the published observations.
Five external software/archive/binary pins remain provenance metadata in this
mode; their bytes are not claimed to have been independently rechecked.
It requires no CNS installation or network access and writes no scientific
artifact. Full local release verification, when the separately downloaded
archives and installed binaries are available, remains
`python -m aav_benchmark.haddock_method`. The historical
`audit_reproduction.py` generates the original report using those local tools
and is not the portable replay command.

Fourteen tests passed. They check actual source pins against the pinned release
archive, atom rebuilding, physical topology/neutrality, glycosidic links,
handedness, AIR information, and 20→12→12 workflow inventories. Corruption
controls remove atoms, reflect a stereocenter, corrupt a topology bond, remove
required pin fields or relabel the guided run as blind; all are detected. Tests
are read-only and do not execute CNS.

Release-integrity tests require the locally downloaded archives and isolated
installation and explicitly skip when those external tools are absent. Chemical
artifact tests use retained scientific outputs; they do not reinstall software
or launch CNS. The reported fourteen-pass run used the actual installed release.
