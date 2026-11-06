# Archived receptor-neighborhood coverage audit and 8 Å amendment

**Three of 52 frozen docking boxes have omitted source receptor atoms inside
Vina's ordinary 8 Å nonbonded range.** All three are AAV development boxes.
The other nine AAV boxes and all 40 prepared external boxes already cover that
range. Every external box retains its whole archived source receptor, so this
repair changes none of the external atomic-reference inputs.

The [source-only audit](../results/blind_pose/neighborhood_coverage_audit.json)
recomputes coverage from archived receptor inputs and frozen search boxes. It
reads no docking poses, docking scores, pose metrics or accuracy tables. The
[amendment](../benchmarks/blind_pose_neighborhood_amendment.json) preserves the
original 6 Å experiment and freezes a separate three-box correction before
corrected docking. This is development geometry repair and sensitivity work,
not a pristine confirmatory experiment or an atomic-reference qualification.

## Authoritative cutoff and geometric rule

[AutoDock Vina 1.2.7 scoring_function.h](https://raw.githubusercontent.com/ccsb-scripps/AutoDock-Vina/v1.2.7/src/lib/scoring_function.h)
assigns 8.0 Å to the ordinary Vina nonbonded terms and `m_cutoff`. The same source
also has a 20 Å maximum for its special linear-attraction term;
[potentials.h](https://raw.githubusercontent.com/ccsb-scripps/AutoDock-Vina/v1.2.7/src/lib/potentials.h)
restricts that term to macrocycle glue pairs. The independently generated ligand
inputs here contain no G/CG macrocycle glue types. Both complete upstream source
files are archived with SHA256 pins under
`results/blind_pose_neighborhood8/software/`.

The frozen original experiment selected whole protein chains if any source
heavy atom lay within **6 Å Euclidean distance of the closed search box**.
That buffer is insufficient as a general assertion of complete 8 Å Vina
coverage. Whole-chain retention happens to make many individual 6 Å selections
sufficient, which must be checked geometrically rather than assumed.

For source atom coordinate `x`, box center `c` and side lengths `s`, the audit
calculates `||max(abs(x-c)-s/2, 0)||₂`. It groups the distances by chain and
retains every atom of each chain whose minimum distance is at most the declared
cutoff. This is Euclidean distance around faces, edges and corners, rather than
an axis-aligned buffer with square corners. Equality is conservatively included.
The audit records every source chain's minimum distance and an exact nearest
omitted atom witness.

The guarantee concerns receptor interactions for ligand atoms inside the
**declared closed box**, relative to the archived dry protein-only source
receptor. It does not certify biological assembly completeness, a retained-water
model or accessibility outside that source. Returned-pose evaluation against the
full source assembly remains a separate task.

## Complete preparation denominator and source checks

All 27 frozen condition/method preparation manifests were read. Their prepared
jobs agree exactly with the **52 archived box-input manifests**: 12 AAV and
40 external. Failed preparation, failed prediction and no-candidate conditions
retain their original statuses and are not counted as successful preparations.

For each prepared box the audit:

- Restores only the **original** site-input revision from its pinned compressed
  archives, checks all restored bytes and binds its original source SHA to the
  frozen experiment configuration. Historical absolute metadata paths are not
  opened.
- Reconstructs both the 6 Å and 8 Å complete-chain sets from receptor coordinates,
  then checks that the 6 Å set and nearest omitted distance reproduce the
  original preparation.
- Checks every archived preparation-input file's compressed and restored SHA,
  and compares the prepared PDB and PDBQT heavy-atom identity sets against every
  heavy atom of the selected source chains.
- Verifies source coordinates and unchanged independently embedded ligand bytes.
  **All 52 PDB and PDBQT heavy-atom coordinate comparisons are exact: maximum
  displacement 0.0 Å.** No source heavy atom was added, lost or moved.

The geometric algorithm and original configuration hashes are recorded in the
report. No original configuration, preparation, runner, helper or execution
output was modified.

## Three affected boxes

| Frozen box | Nearest omitted atom at 6 Å selection | Added whole chains for 8 Å | Added heavy atoms | Nearest omitted atom after correction |
|---|---|---|---:|---:|
| `aav9_apo/oracle/box01` | E:LYS664:CE, 7.334804 Å | B, E | 8,262 | 17.713127 Å |
| `aav9_apo/p2rank/box03` | Y:PHE662:CD1, 6.963908 Å | Y | 4,131 | 13.884203 Å |
| `aav9_holo/p2rank/box05` | G:PHE662:CD1, 7.484907 Å | G | 4,131 | 14.239322 Å |

All three corrected receptors retain five complete chains and 20,655 source
heavy atoms. The source audit determines which inputs need correction; no
success, failure, energy or recovered-pose property enters that selection.
There is no qualified AAV atomic reference in this experiment, and the coverage
repair does not create one.

The remaining **49/52 inputs** already contain every source chain required by
the 8 Å rule. This is not a reason to relabel the original 6 Å policy as generally
8 Å complete. Both the original policy and each box's verified coverage remain
visible.

## Separate correction and execution boundary

The correction preserves the original independent ligand PDBQT bytes, box
center and dimensions, candidate ranks, protein-only dry rigid model, three
seeds (`104729`, `130363`, `155921`), two CPUs, exhaustiveness 16, 20 returned
modes, 0.5 Å pose-separation setting and 10 kcal/mol output window. Only the
source-chain coverage rule changes from 6 to 8 Å. Full-chain preparation uses
the same audited Meeko routine, with its heavy-atom preservation checks.

New inputs live under `data/processed/blind_pose_neighborhood8/`; their portable
compressed copies, separate configuration, freeze receipts and later raw runs
live under `results/blind_pose_neighborhood8/`. The correction manifest maps
all three boxes to their original input hashes. Original results remain intact;
a later report must identify which variant it evaluates and must not silently
replace legacy observations.

A complete local package and driver snapshot is pinned under the correction's
`software/` folder. The new driver implements 8 Å selection directly and calls
common preparation/execution functions. It does not monkeypatch or edit any
legacy helper. Running through that snapshot preserves the actual preparation
and execution code if repository helpers subsequently change.

The amendment is frozen, all three corrected input archives are prepared, and
[offline preflight](../results/blind_pose_neighborhood8/preflight_receipt.json)
passed. Their PDB and PDBQT heavy coordinates also match source exactly. The five
common helper files are byte-identical to the original execution-start snapshot.
**Nine corrected Vina invocations may start only after
the coordinator grants a CPU slot.** Source auditing, freezing and preparation
commands launch no docking. The new untouched workflow must use 8 Å coverage
before its first freeze; these development corrections cannot be described as
such an untouched test.

## Offline reproduction

Audit the original 52 input archives without running predictors or docking:

```sh
python scripts/audit_neighborhood_coverage.py audit --check
python -m pytest -q tests/test_neighborhood_coverage.py
```

After the correction inputs have been prepared, verify every archived corrected
input against the original source receptor and unchanged ligand bytes:

```sh
python results/blind_pose_neighborhood8/software/scripts/audit_neighborhood_coverage.py \
  --root "$PWD" verify
```

The analytic tests cover interior points, box faces/edges/corners, inclusive
cutoffs, whole-chain retention, source witnesses, translation invariance,
nonfinite boxes, duplicate atom identities and invalid cutoffs. They do not
interpret a software execution as biological validation.
