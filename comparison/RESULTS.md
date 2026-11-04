# AAV9 apo versus holo: archived local-docking comparison

Six completed AutoDock Vina 1.2.7 runs compare ligand-stripped holo (7MTZ)
with independently solved apo (7MT0), using the same three seeds.
This is a retrospective coordinate diagnostic for one biological system.

## Results

Values below are seed means. Displacement is named-heavy-atom RMSD from
one deposited modeled GAL copy, in the aligned receptor frame.

| State | Runs | Retained poses | Top1 score (kcal/mol) | Top1 displacement (Å) | Top5 minimum (Å) | Best retained (Å) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| holo | 3/3 | 60 | -4.866 | 11.494 | 8.470 | 3.333 |
| apo | 3/3 | 60 | -4.617 | 9.055 | 5.098 | 3.770 |

**0/6 runs** retained a pose within 2 Å of the chosen modeled reference.
Favorable Top1 Vina scores coincide with substantial coordinate displacement
in both states. These observations do not establish experimental binding
affinity, specificity, or biologically incorrect recognition.

## Individual runs

| State | Seed | Poses | Top1 score (kcal/mol) | Top1 displacement (Å) | Top5 minimum (Å) | Best retained (Å) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| holo | 104729 | 20 | -4.865 | 11.491 | 8.474 | 3.333 |
| holo | 130363 | 20 | -4.869 | 11.487 | 8.465 | 3.339 |
| holo | 155921 | 20 | -4.863 | 11.505 | 8.471 | 3.328 |
| apo | 104729 | 20 | -4.620 | 9.045 | 5.095 | 3.681 |
| apo | 130363 | 20 | -4.619 | 9.047 | 5.105 | 3.724 |
| apo | 155921 | 20 | -4.611 | 9.074 | 5.096 | 3.905 |

## What was compared

Both inputs retain three whole subunits (author chains A, B and I) around
the same reference-derived 20 Å cube. The historical `complete` suffix
means this local neighborhood; the docking input is not a full capsid.
Before cropping, apo coordinates were aligned using 31,080 corresponding Cα
atoms from the full assemblies (0.607 Å RMSD). The ligand was not fitted separately.

The archived setup used standard Vina scoring, three seeds (104729, 130363,
155921), exhaustiveness 16, two CPUs and up to 20 returned poses per run.
It used a rigid receptor and dry preparation. Input ligand coordinates were
generated independently; the search box and neighborhood used the known site.

## Interpretation and limits

The apo input has a lower mean Top1 coordinate displacement in these runs;
the holo input has a more favorable mean Top1 Vina score. Neither quantity
establishes which state models experimental recognition better.

The AAV9 modeled GAL coordinates have not passed an independent atomic-pose
reference qualification here. RMSD and the 2 Å count are descriptive coordinate
agreement, not validated pose accuracy. This local, known-site experiment
does not measure blind site discovery. Seeds are computational repeats, not
independent biological samples. One fixed reference copy is used; no capsid-wide
symmetry minimization or internal sampling diagnosis is claimed. There are
no tested nonbinders or assay measurements in this comparison.

The wider hypotheses in [PROTOCOL.md](../PROTOCOL.md) remain research proposals.
This completed comparison answers the smaller descriptive apo/holo question.

## Reproduction and sources

Run `python3 pilot/pilot-script.py --check` from the repository root.
It verifies input hashes and recomputes both this table and `results.json`
from the archived PDBQT models, using Python's standard library. It does not
rerun docking or independently reproduce receptor/ligand preparation.

[Input provenance](inputs/SOURCES.json) records each original path, SHA256
and extraction. Only the six selected runs and the reference fields needed
for this calculation were extracted; poses and logs are unchanged.
Source snapshot: `ddd84f9c9021d43a4c58868c9116d0523f3cbb79`.
