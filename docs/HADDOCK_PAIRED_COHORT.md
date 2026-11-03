# Prepared paired HADDOCK feasibility workflow

`aav_benchmark.haddock_paired_cohort` prepares the concrete source-only input plan for six fixed parents, two protein preparations, three site-information arms and three seeds. It has not generated a cohort docking pose. An external coordinator freeze matching the exact prepared plan is required by both native-execution entrypoints.

The new protocol is `topoaa → rigidbody(20) → flexref(all available initial models, at most20)`. All final models are retained and ranked by individual HADDOCK score. There is no clustering, interface residue selection, native reference or internal pose-accuracy computation. This explicitly adapts the released glycan software test: its actual `20→12` retention arose from native-informed clustering and was not a published fixed12 parameter. The new protocol retains `0..20` refinements and makes no claim to reproduce the original clustering workflow or recommended full sampling. The released [glycan test configuration](https://raw.githubusercontent.com/haddocking/haddock3/815fcc14100896e35d1795058465cb547f37a8c5/examples/docking-protein-glycan/docking-protein-glycan-test.cfg) remains the source of the20-model test scope, rigidbody `w_vdw=1` and tolerance20.

## Frozen inputs and explicit denominator

`results/structural_pose_v3/vina_carb_prepared_inputs.json` (SHA256 `0dc9e3975d5451d3e21b2f474b7218ecf4d82e110e197fab95630bb13c22579b`) supplies all36 condition identities and pins their source-only `condition.json` files. The adapter reads the original rank-one candidate from `candidate_selection` for fpocket/p2rank, or the already declared oracle center from its first job. It takes box sizes from the corresponding v3 independent parent chemistry. It never reads a Vina pose, run result, recovery value or VinaCarb score. A Vina-specific preparation failure therefore does not erase a valid HADDOCK input: all six4R9F site/profile conditions are available because actual CNS supports its unchanged MSE-containing protein.

All108 seed-level cells remain in `results/haddock_paired_cohort/prepared_plan.json`:

| Source-only status | Conditions | Seed cells |
| --- | ---: | ---: |
| Ready for native execution | 13 | 39 |
| Exact protein-input parameterization gate failed | 15 | 45 |
| Bounded adapter does not preserve multi-chain receptor assembly (3ZWE) | 6 | 18 |
| Actual p2rank output has no candidate (2RDK, both protein profiles) | 2 | 6 |
| Total retained | 36 | 108 |

The complete source-only topology audit and its conservative failure matrix are documented in `HADDOCK_COHORT_PREPARATION.md`. Its immutable preparation receipt is SHA256 `9acaf59029d153b66764bacb8ba485ad26b7c27cc097b57237e1e2d8f62756cf`. Initial exact receptor topology support is strict2RDK/3ZWE/4R9F and repair2RDK/3ZWE/4R9F/5V6F. A separate pre-pose assembly-protocol gate excludes3ZWE from this bounded adapter; the final native-ready profiles are strict2RDK/4R9F and repair2RDK/4R9F/5V6F. CNS heavy-atom rebuilding or its conditional protein minimization prevents an exact matched-input claim for the other five case/profile combinations; these failures are retained, without native-informed correction or substituted method inputs.

All prepared inputs retain the same complete v3 independently initialized parent, every observed protein assembly chain and literal source protein chemistry. CNS glycan names are standard CCD names and the ligand segid is `L`; `data/preflight_preparation.json` maps every named heavy atom to its canonical `branch:CCD:atom` key. The independent ligand's rigid parking translation prevents accidental topology-stage covalent glycosylation and leaves its intraligand geometry unchanged. Observed protein maps remain distinct from modeled repair atoms. Post-docking protein flexibility is method behavior; it must not be confused with the source-input equality gate.

## Pre-pose multi-chain protocol amendment

`benchmarks/haddock_assembly_protocol_amendment.json` records an explicit source-only amendment made before any HADDOCK cohort pose. Initial topology preserved the3ZWE trimer, but stock [input grouping](https://raw.githubusercontent.com/haddocking/haddock3/815fcc14100896e35d1795058465cb547f37a8c5/src/haddock/libs/libcns.py) counts unique chain/segids as components. [Rigidbody](https://raw.githubusercontent.com/haddocking/haddock3/815fcc14100896e35d1795058465cb547f37a8c5/src/haddock/modules/sampling/rigidbody/cns/rigidbody.cns) separates and moves each independently, so A/B/C+L would be four bodies. Default [flexref automatic selection](https://raw.githubusercontent.com/haddocking/haddock3/815fcc14100896e35d1795058465cb547f37a8c5/src/haddock/modules/refinement/flexref/cns/flexauto-neighbors.cns) also flexibilizes protein–protein interfaces, with [rigid groups per segment](https://raw.githubusercontent.com/haddocking/haddock3/815fcc14100896e35d1795058465cb547f37a8c5/src/haddock/modules/refinement/flexref/cns/torsiontop.cns). Rigidbody origin flags alone do not freeze flexref atoms, and nseg=0 does not join separate chains.

The uniform operational rule therefore requires one receptor chain/segid for this bounded native-default adapter. It does not remove a chain, merge source identities, invent native glycan restraints or introduce an unvalidated grouping plugin. All six3ZWE conditions and18 seed cells remain explicit `assembly_protocol_unsupported` failures. This describes an adapter/workflow limitation under the preserved-assembly requirement, **not** a claim that HADDOCK cannot dock multi-chain receptors; its official antibody workflow uses interchain restraints. A validated source-protein-only grouping/restraint adaptation would be a separate future protocol. The Vina/VinaCarb full-trimer study is unchanged.

## Exact box-to-AIR adaptation

The adapter selects exactly the protein heavy atoms inside the closed frozen box: each absolute coordinate difference from its supplied center must be no larger than half that axis's full side length. It verifies each selected `(segid,resid,name)` against the actual CNS protein topology. It writes one OR-union of these atom identities against one OR-union containing every full-parent glycan heavy atom, with `2.0 2.0 0.0` distance parameters. The [official restraint documentation](https://www.bonvinlab.org/haddock3-user-manual/intro_restraints.html) defines this selection and pseudo-distance syntax. Empty selections fail explicitly. There is no residue expansion, contact calculation, surface filter, nearest-atom fallback or box enlargement.

This is an ambiguous patch-contact prior, rather than a hard Cartesian search box. It is not exact information equivalence to Vina's search-space restriction. The two predicted arms use their own unchanged protein-only predictor centers; the oracle arm is explicitly native-derived box-center information. Neither predicted arm gains a native-contact restraint. `box_air.json` retains exact atom identities, box dimensions, source/topology hashes and AIR checksum.

Every one of the13 final prepared AIRs was read by the real pinned CNS binary against the audited receptor and glycan topologies. The check only loads topology/coordinates and prints a restraint summary: no atom is moved, optimized or written. All13 yielded exactly the expected protein/glycan selection counts and one NOE, without a parser error. The largest final union selects830 protein heavy atoms and23 glycan heavy atoms; the largest eligible full parent has34 heavy atoms. Earlier source-only19-AIR syntax evidence included3ZWE, but successful AIR parsing does not establish assembly preservation. This avoids assuming that a large generated OR expression is valid merely because the text looks plausible. Scripts, output and the receipt are under `AIR_syntax_preflight/`.

## Native settings, scores and execution guard

Each job uses one CPU, explicit base seed104729/130363/155921, topoaa seed917 and the source-only observed-OXT terminal rule. Rigidbody uses sampling20, tolerance20 and `w_vdw=1`. Flexref has sampling factor1, tolerance20 and the pinned default MD step counts500/500/1000/1000. Random removal of the single AIR is disabled in both stages. The pinned [rigidbody](https://raw.githubusercontent.com/haddocking/haddock3/815fcc14100896e35d1795058465cb547f37a8c5/src/haddock/modules/sampling/rigidbody/__init__.py) code gives model seeds `base+1..base+20`; the [flexref](https://raw.githubusercontent.com/haddocking/haddock3/815fcc14100896e35d1795058465cb547f37a8c5/src/haddock/modules/refinement/flexref/__init__.py) code inherits model seeds at sampling factor1. Failed slots are not replaced and the pool is never enlarged.

The full108-cell nominal ceiling is2160 initial and2160 refined models; current source/protocol gates allow at most780 of each to be attempted. The plan records actual stage file counts, unfilled initial slots, model atom counts, native seeds, exact file hashes, method energy components and individual full-precision score orders. The precision is that of native HADDOCK model metadata, without additional rounding. Ties use the immutable original model number. Every physical model is retained, including unranked/error outputs; incomplete final score metadata makes the execution fail. Cluster fields stay null because no clustering takes place. These model-file/score records do not calculate experimental accuracy.

The plan pins the official distribution, unchanged CNS binary, release lock, native environment freeze, all370 installed Python/CNS/topology/parameter/default source files, this helper, proposal, source manifests, each input, box/AIR and config. Three real `haddock3 --setup` checks passed without native execution: the3ZWE trimer/charged terminus,4R9F MSE and repaired5V6F cases. Eight regression tests cover the full denominator, method-specific source eligibility, actual CNS AIR parsing, strict closed-box selection, no residue fallback, source-center reuse through Vina-specific failures, bounded configs, external freeze/pin rejection and complete stable individual score export.

Root must provide an external JSON receipt with `status="frozen_before_first_HADDOCK_cohort_pose"`, `execution_authorized=true`, and `prepared_plan.sha256` matching the final plan. The helper additionally checks its own source and all upstream/runtime pins. It refuses an existing or partial native run; it does not silently restart or replace models. The foreground traversal writes terminal progress after each of all108 cells. Completed matching job receipts may be reused on resumption; a partial native directory requires explicit coordinator adjudication.

```sh
# Source-only preparation into a fresh output tree:
.venv/bin/python -m aav_benchmark.haddock_paired_cohort --prepare-only

# Coordinator only, after the separate freeze, sequential foreground/one CPU:
.venv/bin/python -m aav_benchmark.haddock_paired_cohort \
  --run-all-foreground \
  --plan results/haddock_paired_cohort/prepared_plan.json \
  --freeze PATH_TO_ROOT_FREEZE.json
```

The public preflight evidence can be restored before replaying tests using:

```sh
tar -xzf results/haddock_prospective_preparation/topology_artifacts.tar.gz \
  -C results/haddock_prospective_preparation
```

The new source-only proposal is `benchmarks/haddock_paired_cohort_v2_proposal.json`; the earlier proposal remains archived with its preflight history. A separate external metric adapter must map CNS CCD identities through the canonical keys to the existing qualified multi-copy atom masks. Experimental native-pose metrics are permitted only after every planned cell is terminal and the independent external analysis code/plan has been frozen. No cohort accuracy has been read or computed by this helper.
