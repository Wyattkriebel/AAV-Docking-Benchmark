# P2Rank method and reproducibility

This adapter executes the official **P2Rank 2.5.1** release, using its shipped
`default` random-forest prediction model. The version, model, configuration,
runtime, and operational limits were fixed before development-cohort predictions.
The official [release page](https://github.com/rdk/p2rank/releases/tag/2.5.1)
identifies 2.5.1 as the latest stable release when checked on 2026-10-30;
2.6-alpha was a prerelease. The release states that its models are unchanged
from 2.4. This is P2Rank prediction, not rescoring of another method's sites.

`benchmarks/p2rank.lock.json` pins the official archive URL and SHA256, source
commit `9808a7723be9a94e2ffc21ab5f724cb6ae4ba01e`, MIT license, 136 production
resource checksums, and a portable Eclipse Adoptium Temurin **21.0.13+9-LTS**
macOS aarch64 JRE. The JRE is GPLv2 with the Classpath exception; its archive and
the Java executable/JVM library are pinned separately. No global Java or
package-manager installation is required. P2Rank 2.5.1's release notes support
Java 17–24. Other platforms/runtimes require a reviewed lock update.

## Frozen execution

The adapter uses this command, with absolute input/output paths:

```text
prank predict -f receptor.cif -o raw -threads 1 -visualizations 0 -seed 42 -model default
```

All other prediction parameters retain the pinned release defaults. Disabling
visualizations avoids copies and rendering artifacts; it does not tune pocket
prediction. The unmodified launcher limits the Java heap to 2048 MiB. Each
receptor has a 600-second prediction timeout; the adapter terminates the Java
process group on timeout. Version/integrity checks and parsing also contribute
to the recorded total wall time. Receptors run individually with one worker.
User-provided JVM option variables are removed to prevent hidden changes.

The shipped default feature set includes B factors. Upstream recommends the
`alphafold` profile for AlphaFold, NMR, and cryo-EM inputs, but this benchmark
fixes one default model across receptors. AAV cryo-EM performance therefore has
that declared limitation; the profile is not changed after inspecting outcomes.
The model predicts small-molecule ligand sites. Exposed glycan-binding patches
and a full viral capsid are not established in-domain validation settings.

## Input and output contract

```python
from aav_benchmark.p2rank import run_p2rank

result = run_p2rank(
    "full_reference_free_receptor.cif", "fresh_output_directory",
    "tools/p2rank_2.5.1/prank",
)
```

Both PDB and mmCIF are accepted. The caller supplies the complete prepared
protein receptor/assembly. The adapter reads only that receptor and the frozen
method installation. It does not read a native ligand, reference center,
contact list, docking box, or ligand-derived crop. For the full AAV assembly,
mmCIF avoids PDB atom-serial overflow. Numeric mmCIF entity IDs are necessary
for the bundled BioJava reader; long author chain IDs are supported.

Each pocket contains the upstream rank, raw score, calibrated probability,
Cartesian center in Å, SAS-point count, upstream surface-atom IDs, and residue
footprint. Footprints come from `*_predictions.csv`; they may overlap between
pockets. `*_residues.csv` validates that every supplied amino-acid residue was
reported. Residues map exactly to the input's author chain, signed residue
number, insertion code, and residue name. Missing, duplicated, or mismatched
identities fail explicitly. Original label identifiers and assembly-source
mapping remain in the receptor-preparation ledger; upstream surface-atom IDs
are retained as raw method metadata, not asserted to be deposited atom IDs.

`result.json` records the input SHA256, executable and lock hashes, exact
command, Java home, parameters, status, wall time, return code, predictions,
and checksums/lengths for all raw artifacts. Standard output, standard error,
version output, upstream CSVs, `params.txt`, and upstream logs are retained.
An empty prediction table after a successful complete run is `status: ok`
with zero pockets. Input, installation, execution, output-validation, and
timeout failures retain an explicit result and their runtime; they remain in
the unconditional receptor denominator. Output directories must be fresh to
prevent accidental reuse of stale predictions.

## Model provenance and overlap audit

The release's [model inventory](https://github.com/rdk/p2rank/blob/9808a7723be9a94e2ffc21ab5f724cb6ae4ba01e/distro/models/readme.md)
states that the default model was trained on `chen11-fpocket.ds` using
`default.groovy`. The [dataset repository](https://github.com/rdk/p2rank-datasets/tree/0236ecb38cbb60b89849a1ea36fbd9ee93f3e906)
lists training and evaluation collections. This audit uses the **dataset
definition files**, not every structure stored in their directories, as the
official dataset README requires.

At dataset commit `0236ecb38cbb60b89849a1ea36fbd9ee93f3e906`, case-insensitive
exact PDB-ID checks found no matches for 154L, 1GNY, 1KJL, 1QFO, 2J1T, 1PWB,
2VUZ, 7MTZ, or 7MT0 in these lists:

| Dataset | Entries | Unique PDB IDs | Exact matches |
|---|---:|---:|---:|
| CHEN11 / CHEN11-fpocket | 251 each | 241 each | 0 |
| JOINED | 537 | 537 | 0 |
| HOLO4K | 4009 | 4009 | 0 |
| FPTRAIN | 222 | 222 | 0 |

The lock stores the exact URLs and content checksums for these audited lists.
This is an exact structure-identifier audit. Sequence/family overlap has not
been excluded. The inspected provenance also does not establish independent
feature-table or score-calibration data. The absence of exact PDB matches does
not justify a held-out or family-generalization claim.

## Installation and validation

On the pinned macOS arm64 platform, installation can be replayed with
`aav_benchmark.p2rank.install_p2rank()`. It downloads and checks the frozen
archives, extracts into local `tools/`, and verifies the installation and runtime.
Existing tool directories are verified without replacing their contents.

Download the two archive URLs in the lock into `tools/` and verify their
SHA256 values before extraction. Extract P2Rank as `tools/p2rank_2.5.1` and
the JRE as `tools/jre-temurin-21.0.13+9` (the archive's original top directory
is `jdk-21.0.13+9-jre`). Preserve the release files without modifying the
launcher or configuration. Archives and installations are local ignored
tools, not vendored repository artifacts. Compressed downloads total about
324 MB and expanded files about 462 MB.

Eight adapter tests cover original residue identities, zero predictions,
overlapping footprints, incomplete/mutated maps, nonfinite output, invalid
ranks, stale-output rejection, archived installation failure, and process-group
timeout cancellation. Before cohort execution, actual official P2Rank runs
also completed on a 100-ATOM ligand-free fragment in PDB and mmCIF formats,
including a two-character author chain ID. Both completed with zero pockets;
these are runtime/parser smoke checks, not site-prediction performance evidence.
No development receptor was run during those smoke checks.
