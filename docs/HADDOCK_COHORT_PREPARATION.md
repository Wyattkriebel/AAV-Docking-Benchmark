# Source-only HADDOCK cohort preparation and unfrozen pilot proposal

All six new independent v3 full glycans have passed actual HADDOCK3/CNS parameterization. This establishes chemical-input feasibility beyond the earlier 1LMQ released-example software test. It does not establish cohort docking performance. No `rigidbody`, `flexref`, AIR-guided pose generation, pose scoring or accuracy evaluation has been run for this work.

The unchanged release remains HADDOCK3 `2026.10.1`, source commit `815fcc14100896e35d1795058465cb547f37a8c5`, with the official arm64 CNS binary SHA256 `e6febf74eccb740f8bc5d49955cc8fc202025d5d642e129e7738e7280ad95a88`. The original `benchmarks/haddock.lock.json`, released-example results and helper remain unchanged. Its documented licensing context still applies; this work performs no registration or additional legal affirmation.

## Complete chemical parameterization

`aav_benchmark/haddock_cohort_preparation.py` reads only the frozen v3 independent full-parent SDFs and their source-only reports. It converts canonical keys `branch:CCD:atom` into the same CCD residue and atom names for CNS; no GLYCAM alias assumption, native-bound coordinate, new embedding or free minimization enters the input. Every input heavy atom and full glycosidic branch graph is retained. Input hydrogens are deliberately omitted and explicitly rebuilt by CNS.

The policy `benchmarks/haddock_prospective_policy.json` was written before the six actual CNS jobs. `results/haddock_prospective_preparation/execution_plan.json` pins all six input PDBs, configurations, source driver, policy and CNS binary before the first job. Configurations contain only `[topoaa]`, `ncores=1`, `debug=true`, `clean=false`. Jobs ran sequentially, without receptors or docking stages. The [official topology documentation](https://www.bonvinlab.org/haddock3-user-manual/modules/topology.html) explicitly says topology generation can rebuild missing atoms; residue-name presence in `carbohydrate.top` is therefore insufficient evidence of support.

The audit requires exact named heavy inventories, actual CNS heavy-bond connectivity including every interresidue link, real atom types/masses, total partial charge consistent with the full parent's formal charge, no missing/rebuilt heavy atom, coordinate changes ≤0.002 Å, every v3 ring/named-volume gate, and a fresh complete isomeric-graph assignment from the post-CNS 3D coordinates. Six actual outputs pass:

| Parent | Input/output heavy atoms | CNS-built hydrogens | Actual heavy bonds | Maximum heavy movement |
| --- | ---: | ---: | ---: | ---: |
| 2RDK | 23 / 23 | 22 | 24 | 0 Å |
| 3ZWE | 44 / 44 | 42 | 47 | 0 Å |
| 4R9F | 23 / 23 | 22 | 24 | 0 Å |
| 5HZA | 33 / 33 | 32 | 35 | 0 Å |
| 5V6F | 34 / 34 | 32 | 36 | 0 Å |
| 8I5U | 23 / 23 | 22 | 24 | 0 Å |

All are neutral to floating-point summation precision; all six CNS coordinate-read logs contain zero `READC-ERR` messages. Standard names alone did not grant acceptance: actual PSF connectivity, charges, output coordinates and full stereochemistry were inspected. Each case retains its preparation, exact config, CNS input, log, prepared PDB, PSF and `topology_audit.json`. Six regression tests replay all actual outputs and reject lost/moved heavy atoms, removed heavy bonds, incorrect total charge, absent chemical types and source-byte tampering.

## Paired receptor preflight and limits

`scripts/preflight_haddock_pairs.py` adds the exact strict or repair V3 receptor PDB from `data/processed/heldout_pose_v3/{strict,repair}/CASE`. The receptor bytes and every assembly chain are retained, including the 3ZWE trimer; observed protein atom maps and modeled repair maps remain separate. Glycan segid `L` is reserved and collision-checked. Before topology only, the independent glycan is rigidly translated so its minimum x coordinate lies 50 Å beyond the receptor's maximum x; all glycan hydrogens are translated consistently in a provenance SDF. This prevents CNS distance-based glycosylation detection from introducing an accidental protein–glycan covalent bond. It preserves all intraligand geometry and uses no native glycan geometry.

All 12 strict/repair topology jobs completed with passing complete-glycan gates. Protein output is audited separately against its exact input: every heavy name, every chain, shared coordinate changes, missing/built atoms, CA signed-volume preservation and actual CNS total partial charge are recorded. Protein charge is the method's explicit protonation/topology outcome; a complete formal-charge equivalence to the Vina preparation is not asserted.

The released default `charged_cter=false` removes observed terminal OXT atoms in 3ZWE and 5V6F. The explicit source-only amendment `benchmarks/haddock_terminal_amendment.json` applies a uniform chemical rule before any cohort docking: if every source chain has one OXT at its last observed residue, use `charged_cter=true`; if no source chain has OXT, retain false; mixed, duplicate or nonterminal inventories fail. The four changed case/profile combinations are reparameterized; unchanged configurations reuse their v1 executions, and all 12 remain in the amendment summary. Original default results are preserved.

The default source code also conditionally performs a 100-step Powell minimization when bonded energy exceeds 10,000: [pinned `generate-topology.cns`](https://raw.githubusercontent.com/haddocking/haddock3/815fcc14100896e35d1795058465cb547f37a8c5/src/haddock/modules/topology/topoaa/cns/generate-topology.cns). This changes observed protein coordinates in 5HZA and 8I5U even after the independent root repair. Such inputs cannot be described as the same fixed receptor geometry used by Vina. The strict arm additionally builds 39/7/13 missing heavy atoms for 5HZA/5V6F/8I5U. These are real preparation differences, not software-installation failures or evidence of biological nonbinding. No source-only tuning loop is used to hide them. The complete case/profile denominator is retained; exact matched-input availability and default HADDOCK chemical support are distinct fields.

The actual topology evidence is under `results/haddock_prospective_preparation/paired_receptors/`, with the observed-terminal amendment under `observed_termini_v2/`. A future matched-input study must retain unavailable source-gate cells explicitly. A separate study could intentionally permit HADDOCK's own receptor remodeling, but would need a new declared preparation condition and could not claim identical receptor geometry.

After the observed-terminal amendment, the final source-only availability is:

| Parent | Strict exact receptor | Repair exact receptor | Remaining difference |
| --- | --- | --- | --- |
| 2RDK | Yes | Yes | None |
| 3ZWE | Yes | Yes | Three observed OXT atoms retained by the declared charged-terminus amendment |
| 4R9F | Yes | Yes | None |
| 5HZA | No | No | Strict builds 39 heavy atoms and moves observed atoms ≤0.984 Å; repair moves observed atoms ≤2.647 Å |
| 5V6F | No | Yes | Strict builds seven heavy atoms; observed OXT is preserved by the amendment |
| 8I5U | No | No | Strict builds 13 heavy atoms and moves observed atoms ≤0.996 Å; repair moves observed atoms ≤0.997 Å |

Thus seven of 12 case/profile inputs satisfy exact receptor geometry, with all six case identities and all 12 source profiles retained. All 22 actual CNS jobs (six glycans, 12 default pairs, four terminal-amended pairs) returned zero. Twelve source/actual-output/corruption tests pass. The model-free source artifact receipt is `results/haddock_prospective_preparation/preparation_receipt.json`.

## Proposed information-matched feasibility pilot

`benchmarks/haddock_paired_cohort_proposal.json` is **unfrozen** and authorizes no docking. It retains six parents, strict/repair source profiles, frozen predicted-top1/oracle box arms and three base seeds: 72 workflow cells before source-gate exclusions. All inputs come from the same v3 independent full parents and frozen V2/V3 protein/box sources. Native-contact AIRs and native-interface clustering lists are prohibited.

The proposed AIR generator selects exactly the protein heavy atoms inside the closed frozen cube, using its supplied center and full side lengths, then verifies their actual CNS `(segid,resid,name)` correspondence. It writes a single OR-union of those atoms against all independently prepared glycan heavy atoms: `assign (protein_union) (glycan_union) 2.0 2.0 0.0`. No residue expansion, solvent-accessibility filter, nearest-atom fallback or empty-box enlargement is permitted; empty selections fail. This avoids asserting that every box residue must contact the ligand. Selected atom identities and source/box/AIR hashes must be frozen. The [official restraint syntax](https://www.bonvinlab.org/haddock3-user-manual/intro_restraints.html) describes the pseudo-distance and selection syntax; the [pinned restraint implementation](https://raw.githubusercontent.com/haddocking/haddock3/815fcc14100896e35d1795058465cb547f37a8c5/src/haddock/modules/sampling/rigidbody/cns/read_noes.cns) uses ambiguous averaging. This is a coarse spatial prior rather than a hard Cartesian box, so it must not be called exactly equivalent to Vina's search restriction.

The concrete proposed stages are `topoaa → rigidbody(20) → flexref(all surviving models, sampling_factor=1)`. Retain released test settings `tolerance=20`, rigid-body `w_vdw=1`, and default flexref step counts `500,500,1000,1000`. Disable random AIR removal in both stages so the single site restraint cannot disappear. Omit every native-reference `caprieval`, native-informed RMSD clustering selection and native-interface residue list. Refining all 20 instead of applying those released clustering/selection stages is an explicit prospective adaptation, not a claim to reproduce the exact released workflow.

Three base seeds `104729,130363,155921` yield model seeds `base+1` through `base+20`; `flexref` with sampling factor one inherits each model seed, as implemented in the pinned [rigidbody](https://raw.githubusercontent.com/haddocking/haddock3/815fcc14100896e35d1795058465cb547f37a8c5/src/haddock/modules/sampling/rigidbody/__init__.py) and [flexref](https://raw.githubusercontent.com/haddocking/haddock3/815fcc14100896e35d1795058465cb547f37a8c5/src/haddock/modules/refinement/flexref/__init__.py) code. Failed model slots stay failed without replacement. The resulting maximum is 60 refined models per supported case/profile/information arm. This is small-sampling feasibility: the [official full glycan workflow](https://www.bonvinlab.org/haddock3/pages/examples.html#docking-protein-glycan) uses substantially more sampling than its 20-model software test.

Root must freeze actual source receipts, matched-input eligibility, box-derived AIRs, exact configs, seeds, model identities, resources and downstream analysis before generating any cohort pose. The previous 1LMQ example remains software validation using its published native-informed restraints; it does not answer this new cohort comparison.
