# Independent chemical-ideal pyranose initialization

The executable `aav_benchmark.ideal_ring_conformers` produced complete, independently initialized ligands for 1GNY and the six new structural-test parents, without reading experimental ligand coordinates, receptors, docking poses or accuracy outcomes. Every selected SDF passes the same declared chemical and physical gates after serialization. This is a separate input-geometry sensitivity profile; the frozen original preparations and their evidence remain unchanged. It does not authorize a docking study by itself.

The motivation is a source-only physical-input diagnosis: ordinary independent ETKDG/MMFF preparation can leave pyranose rings in boats, twist-boats or an opposite chair. Vina treats these ring coordinates as rigid. A failed downstream glycosidic-geometry term therefore cannot be attributed simply to a component-name support list. The independent Cremer–Pople diagnostic is retained separately in `results/common_pose_rescoring/geometry_diagnosis/`; none of its docking scores or pose accuracy fields was read for this work.

## Chemical inputs and independence

`case_schema()` reads chemical metadata from the unchanged `benchmarks/external_cohort.json` for 1GNY and `benchmarks/heldout_candidate_inventory.json` for the new six. `build_parent_graph()` builds each complete heavy-atom parent from the pinned CCD atom, formal-charge, bond and stereochemical definitions and the declared branch links, removing only explicit leaving atoms. It assigns canonical keys `branch:CCD:atom`, handles the change in tetrahedral neighbor order at substitutions, adds explicit hydrogens, and discards **every monomer conformer before embedding**. The seven complete isomeric graphs match independently extracted, previously frozen preparation chemistry in `tests/fixtures/ideal_ring_expected_parent_graphs.json`.

The only source coordinates used are CCD `_chem_comp_atom.pdbx_model_Cartn_{x,y,z}_ideal` fields for the chemical monomers. The experimental CCD `model_Cartn_*` fields are ignored, as are all protein-bound ligand coordinates. These ideal coordinates are algorithmically generated chemical definitions, rather than an experimental bound conformation; the [wwPDB CCD documentation](https://www.wwpdb.org/data/ccd) describes the ideal-coordinate generation and chemical/stereochemical content. The pinned CCD payload checksums are verified before parsing; every case report records compressed and uncompressed hashes. This work reuses the archived CCD sources and their existing source/license records and downloads no new coordinate source.

The supported scope is one unambiguous six-membered pyranose per component, in named cyclic order `C1,C2,C3,C4,C5,O5`. Unusual ring sizes, incomplete/ambiguous rings, unassigned stereochemistry and unsupported external covalent parents fail explicitly. This initializer makes no claim about seven-membered rings or enzyme-bound nonchair sugars.

## Fixed v3 procedure

The final source-only feasibility policy is `benchmarks/ideal_ring_conformer_policy.json`, profile `ccd_ideal_pyranose_v3`, SHA256 `6d132fc05dc6d391db3457a7b84a788195dae9fac38dc463db96315c194581f4`. The exact executed algorithm is SHA256 `8e1e75f0e67b6ae8376a8b70d3ba02ff4f5d9f5c327e4b9b1f4180641c48321e`; both files and the shared ring-math source are copied into the v3 artifact directory.

Each parent receives exactly 20 candidate slots, with seeds `20261024 + slot`, one CPU and RDKit ETKDGv3, enforced chirality and small-ring torsions. Each solver call has a fixed ceiling of 1,000 internal embedding trials. Those internal retries are distinct from the 20 candidate conformers. There is no adaptive candidate-pool enlargement, stopping after a favorable candidate, native-coordinate map or positional restraint. The [RDKit Book](https://www.rdkit.org/docs/RDKit_Book.html) documents ETKDG and its chirality/coordinate-map options.

After ordinary independent embedding, deterministic MMFF94 minimization restrains all 15 distances within each individual six-atom ring to the CCD ideal distances ±0.025 Å (force constant 1,000), and its six consecutive signed ring torsions to the CCD ideal values ±5° (force constant 100). The maximum is 2,000 MMFF iterations, force tolerance 0.0001 and energy tolerance 0.000001. All constraints are intraring, and each ring can translate and rotate independently; no restraint joins different rings or fixes their position. Absolute signed torsions preserve the handed chemical geometry that pair distances alone cannot distinguish. The [RDKit force-field API](https://www.rdkit.org/docs/source/rdkit.ForceField.rdForceField.html) defines these absolute-distance/dihedral constraints and the minimizer convergence return code.

Every candidate must converge and preserve the complete named heavy-atom graph, formal charges and full-parent isomeric stereochemistry. Every CCD-declared tetrahedral center must retain its named signed volume with the correct sign, absolute magnitude ≥0.05 Å³ and magnitude ≥20% of the corresponding CCD ideal volume. This rejects reflected and nearly planar centers even if an atom's stored chiral tag is unchanged.

Every ring must satisfy all of these operational acceptance rules:

| Measurement | Maximum permitted deviation or range |
| --- | --- |
| Full named Cremer–Pople vector angle from CCD ideal | 15° |
| Named θ difference from CCD ideal | 15° |
| Distance from nearest chair pole | 20° |
| Puckering amplitude Q difference | 0.10 Å |
| Ring bond / ring angle difference from CCD ideal | 0.04 Å / 5° |
| Any intraring pair-distance / signed-torsion difference | 0.08 Å / 10° |
| Absolute ring bond lengths / ring angles | 1.30–1.65 Å / 95–125° |
| Nonbonded heavy separation, excluding graph distances ≤2 | ≥0.65 × sum of van der Waals radii |

These are deliberately conservative development gates for this profile, not universal sugar-conformation standards. A nearest-pole criterion alone accepts opposite chairs; the full **named** puckering-vector criterion prevents that error. `aav_benchmark/ring_geometry.py` implements the original [Cremer–Pople equations](https://doi.org/10.1021/ja00839a011) and separately validates analytical chairs/boats, rigid transforms and the printed paper example. The complete graph and named 3D stereochemistry are also checked after the selected SDF is read back.

All 20 slots are evaluated. Among converged candidates satisfying every gate, selection uses the lowest ordinary MMFF94 energy evaluated at the constrained geometry, breaking ties by slot number. This energy excludes the added restraint energy. It is a deterministic selection criterion within this constrained initialization procedure, not a free-energy calculation or equilibrium conformer population. There is no subsequent free MMFF relaxation and no choice based on experimental similarity, docking or desired outcomes.

## Actual output and preserved development failures

Two prior source-only versions failed before the final v3 search. Both exact policies, complete algorithm sources, all 20 failed slots and logs are retained:

| Version | 1GNY result | Evidence directory |
| --- | --- | --- |
| v1: tight intraring embedding bounds, ordinary initial coordinates | 20/20 failed embedding; 17 initial-coordinate, 3 final-chiral-bounds failures | `results/ideal_ring_conformer_feasibility/profile_v1/` |
| v2: same embedding bounds, random initial coordinates | 20/20 failed embedding; 19 final-chiral-bounds, 1 first-minimization failure | `results/ideal_ring_conformer_feasibility/ccd_ideal_pyranose_v2/` |
| v3: ordinary ETKDG embedding, ring restraints during MMFF | All seven parents have a passing selected conformer | `results/ideal_ring_conformer_feasibility/ccd_ideal_pyranose_v3/` |

The v3 amendment removed the custom embedding bounds, fixed the internal solver ceiling at 1,000, and strengthened graph/nondegeneracy validation before running v3. No final physical gate was relaxed after seeing a v3 result. The version history is source-only algorithm development; it should not be represented as an untouched prospective primary experiment. Policy/report hashes prove the versions executed, and the parent/root messages establish the order; the records do not claim an external preregistration timestamp. No docking or atomic-accuracy outcome was consulted during any version.

| Parent | Full heavy atoms | Rings / stereo centers | Valid / 20 | Selected slot | Maximum named pucker deviation |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1GNY | 46 | 5 / 20 | 2 | 10 | 2.55° |
| 2RDK | 23 | 2 / 10 | 19 | 8 | 3.48° |
| 3ZWE | 44 | 4 / 20 | 2 | 4 | 5.08° |
| 4R9F | 23 | 2 / 10 | 14 | 0 | 0.93° |
| 5HZA | 33 | 3 / 15 | 17 | 14 | 3.53° |
| 5V6F | 34 | 3 / 15 | 13 | 18 | 1.84° |
| 8I5U | 23 | 2 / 10 | 19 | 3 | 4.99° |

The denominator is all 140 fixed slots: 86 valid, 35 embedding failures, 18 stereochemistry failures and one nonconverged MMFF candidate. All failures remain in case reports. Full parents, including atoms outside an atomic-reference scoring mask, are preserved. The 21 selected rings and all 100 selected stereocenters pass. An independent agent re-read all seven saved SDFs, reassigned 3D stereochemistry afresh, compared the frozen complete chemical graphs, and independently checked all gates; its receipt is `results/common_pose_rescoring/geometry_diagnosis/initializer_v3_all7_independent_review.json`.

A complete repeated 2RDK search, using the same 20 seeds without adding candidates, produced byte-identical `report.json` and `ligand.sdf`; it is retained in `ccd_ideal_pyranose_v3/determinism_replay/`. The report hash is `39952bed3406af12e09155f0a6c2e5c8190424bbfd4a329902604eb4303c5146`; SDF hash is `0441e3b260744e50ad10059b82b546fbcfd9bf486e82d8f30d2db5ef9e9b700b`. This establishes reproducibility in the recorded environment, not bitwise identity across RDKit versions or platforms.

## Replay and downstream interface

Generate into a fresh directory; the command refuses to overwrite an existing case report:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  .venv/bin/python -m aav_benchmark.ideal_ring_conformers \
  --output /tmp/independent-ideal-rings --cases 1GNY 2RDK 3ZWE 4R9F 5HZA 5V6F 8I5U
.venv/bin/python -m pytest -q tests/test_ideal_ring_conformers.py
```

The 11 regression tests check all seven independent chemical graphs; pinned-source-only reads; irrelevance of poisoned experimental CCD model coordinates; CCD ideal stereochemical annotations; serialized full-parent geometry; graph/charge/bond/chiral-tag/atom-loss mutations; reflected and planar stereo controls; opposite-chair/boat controls; fixed slots and discarded supplied coordinates; exclusively intraring constraints; and a glycosidic rotation that changes an interresidue torsion while preserving every intraring distance and torsion. The saved-output test also independently reassigns stereochemistry from serialized 3D and verifies minimum-valid-energy selection, source hashes and output hashes.

For downstream preparation, read `CASE/ligand.sdf` with `removeHs=False`. Heavy-atom `canonical_atom_key` properties survive serialization; `CASE/report.json` also supplies `atom_map`, including the RDKit index, branch, CCD name, element and formal charge. Match existing `Lxxx` names by these keys. Preserve this selected geometry through conversion rather than embedding or freely minimizing it again. The report retains every candidate, exact source and policy hashes, templates, complete atom mapping and post-serialization diagnostics. `profile_receipt.json` records the full seven-parent artifact receipt and development history. A separate root freeze must define any downstream study before its first new pose is generated.
