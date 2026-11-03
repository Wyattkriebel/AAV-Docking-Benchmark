# Vina-Carb released-example reproduction

The original authors' **Vina-Carb 1.0, August 30, 2015** distribution now builds
and runs natively on macOS arm64. A fresh installation repeated a released
university tutorial, and separate controls demonstrate that the carbohydrate
intrinsic (CHI) terms are active. This validates this named source/build and
adapter. It does not reproduce the 2016 paper's aggregate performance or qualify
the tutorial coordinates as experimental atomic truth.

The source is the versioned tarball linked from the [official GLYCAM software
page](https://glycam.org/docs/othertoolsservice/downloads/downloads-software/index.html),
with SHA256 `47d170bf434beb74f36d1688658c0ddf715629b79bdbcb0c49da28ed7a68d31f`.
The authors' [installation and input instructions](https://glycam.org/docs/othertoolsservice/download-docs/publication-materials/vina-carb/)
describe the original Vina 1.1.2 base and `VINA_CARB` environment variable.
The current GLYCAM GitHub repository contains later CH–π development; it was
not substituted for the original Vina-Carb release. The distribution carries
an Apache 2.0 license; Boost uses the Boost Software License 1.0. License texts
and exact source/build evidence are retained under
`results/additional_method_reproduction/`.

## Native build and replay

The checksum-pinned [method lock](../benchmarks/additional_method.lock.json)
records the original archive, Boost 1.75.0 source, compatibility replacements,
unchanged scoring-source hashes, input conventions and the auxiliary scorer.
The build uses an isolated Boost installation and C++98, matching the original
code's language generation. It changes neither the CHI functions nor the search
and nonbonded scoring functions.

The compatibility changes replace old Boost filesystem path APIs and provide
the macOS allocator header/name in the bundled GLYLIB header. The allocation
compatibility symbol is absent from the linked executable: the linked GLYLIB
geometry routines do not call it. The original Linux static library is removed
before compiling native objects. The exact replacements and failed initial
build logs remain available. Compiler/library changes can affect native
numerical execution; matching a source checksum is not a promise of identical
binary bytes on every machine.

From the repository's pinned Python environment:

```sh
python -m aav_benchmark.additional_method install --destination /tmp/vina-carb-native
python -m aav_benchmark.additional_method replay-example --installation /tmp/vina-carb-native --destination /tmp/vina-carb-example
python -m aav_benchmark.additional_method build-scorer --installation /tmp/vina-carb-native
```

All installation and replay destinations must be fresh. The installation can
also use `--archives-dir tools` when the two exact source archives are already
cached there. It verifies their hashes before extraction. No system Boost,
shell configuration, installed Java or existing scientific results are changed.
This audited build recipe currently supports macOS arm64. Linux CI verifies
the archived execution and arithmetic without claiming a native Linux build.

## Released tutorial observations

The source example is the [CCBRC protein–glycan docking tutorial](https://github.com/glycodynamics/Tutorial-Protein-glycan-docking/tree/0dc9185a318887a712482650ae667b45c30f1da6),
`flex_lig`, pinned at commit `0dc9185a318887a712482650ae667b45c30f1da6`.
It supplies a prepared humanized Fab receptor, a Lewis Y tetrasaccharide and a
released Vina-Carb output. It uses a known site and a native-like ligand input
conformation. This is a published tutorial replay, rather than a new blind
benchmark case or a complete reproduction of the Vina-Carb paper dataset.

The unmodified configuration uses exhaustiveness 5, eight CPUs, ten modes,
seed 0, CHI coefficient 1 and cutoff 2 kcal/mol. **Explicit seed 0 is fixed in
this original version**: the source calls `auto_seed()` only when the seed
option is absent. An earlier interpretation of that field was corrected by
inspecting the original source; the build-plan snapshot records the original
metadata without changing the executed configuration.

Both native builds returned ten complete poses with 46 heavy atoms per pose
and recognized all three glycosidic links. The second, automated source build
reproduced the first build's pose file byte for byte on this host. Their
execution times were approximately 34.5 and 49.5 seconds; these are observations
under the host's contemporaneous load, not a matched-runtime method comparison.

The released output and the new execution both have top energy −7.6 kcal/mol.
The top-pose coordinate difference between the released computational artifact
and the new artifact is 0.0313 Å. Later rank-matched poses differ by as much as
8.13 Å. The release does not establish the original compiler/runtime, and this
check does not claim equality of all published coordinates. These distances
compare computational artifacts; they are not experimental-reference RMSDs.

## A successful process can still have inactive carbohydrate scoring

The original parser uses atom serials as indices, distinguishes monosaccharides
by residue number, and searches for literal ring names `C1`, `C2`, `C3`, `C4`,
`C5`, and `O5`. Generic `L000` names and one `LIG` residue can therefore produce
a successful Vina-like run without activating the intended carbohydrate terms.

`export_vina_carb_ligand` restores names and sequential residue identities from
the preparation's source atom map. It preserves coordinates, charges, atom
types, serials and the torsion tree. `export_vina_carb_pose` takes coordinates
from a selected retained computational pose while retaining that fixed source
chemical/torsion description. Neither reads experimental ligand coordinates.
The adapter rejects noncontiguous serials, atom loss, unsupported non-C1 donors,
unsupported acceptors and incomplete pyranose rings. It requires the upstream
`VC_log.txt` to confirm every expected linkage; a zero exit status is insufficient.

Atom identity after scoring is restored through invariant serials and the
separate map, because CCD names repeat across sugar residues. Hydrogen identity
is likewise carried by its unchanged serial. Unsupported linkage/ring chemistry
is an explicit method limitation and must remain in the eligible denominator.

A separate probe linked against unchanged scoring objects perturbs one torsion
of the tutorial ligand. On the same coordinates its CHI contribution is:

| CHI coefficient | Contribution, kcal/mol |
| --- | ---: |
| 0 | 0 |
| 1 | 10.4351385876 |
| 2 | 20.8702771752 |

All three links are recognized. This demonstrates active coefficient-dependent
CHI scoring independently of search outcomes or method labels. The probe is
an auxiliary caller; it does not modify the original scoring implementation.

## Scoring a common retained-pose pool

The original `--score_only` branch omits CHI. In addition, it subtracts each
pose's own intramolecular energy, whereas full search uses one common
best-mode intramolecular reference. Adding CHI to `--score_only` would therefore
not reproduce the original within-ligand full-search ranking.

The source-linked auxiliary caller evaluates:

```text
base_shared0(p) = original eval_adjusted(p, intramolecular_reference=0)
combined_shared0(p) = base_shared0(p) + original eval_chi(p, coefficient=1, cutoff=2)
```

The coefficient-0 comparison uses the same `base_shared0` and pose coordinates.
For one fixed ligand and torsion normalization, choosing a different common
intramolecular reference shifts all base scores equally. The returned
`torsion_normalization_scale` must agree before pooling. These energies are
intended for **ranking within one fixed ligand/receptor-context pose pool**;
their shared-zero convention must not be mistaken for a calibrated absolute
affinity or used silently across different ligands.

For all ten retained tutorial poses, this scorer preserves the full-search
rank order. The observed shared-reference offsets span −1.691118 to −1.615010
kcal/mol, admitting one common offset within the full search's 0.1 kcal/mol
printed score unit. This comparison retains coordinate/score rounding rather
than inventing a tighter numerical claim. A separate adapter roundtrip verifies
that restoring source identities leaves all ten scores unchanged to 1e-9.

The ordinary nonbonded cutoff is 8 Å in this source. The five ordinary Vina 1.2.7
terms also use 8 Å in its [pinned source](https://raw.githubusercontent.com/ccsb-scripps/AutoDock-Vina/v1.2.7/src/lib/scoring_function.h).
Its larger maximum includes artificial macrocycle glue terms, which are absent
from these ordinary sugar inputs. Complete neighborhood claims therefore need
an 8 Å interaction audit; a 6 Å margin alone does not prove completeness.

The benchmark comparison requires a separate source and configuration freeze.
This released-example validation does not supply new-cohort accuracy.

## Portable verification

```sh
python -m pytest -q tests/test_additional_method.py
python -m aav_benchmark.additional_method verify-example
```

The latter reads checked-in source inventories, raw tutorial poses/logs and
native scorer output. It checks atom completeness, linkage recognition,
coefficient scaling, rank order and shared-reference arithmetic without a
native executable. Archived raw stdout and its hash accompany every fixed-pose
score. It does not claim to rerun native scoring on Linux or independently
adjudicate experimental density.
