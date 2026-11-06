# Published-site localization results

The eight-system development comparison found the published modeled site among the Top5 distinct candidate classes in 5/8 systems with fpocket and 4/8 with P2Rank. Neither method recovered the AAV9 site in its Top5. These are descriptive localization endpoints; the references qualify approximate published sites rather than validated atomic poses.

## Fixed denominator and primary endpoint

The frozen cohort contains eight systems: AAV9 and seven released glycan–protein complexes. AAV9 apo is a paired ninth condition, reported separately. Three methods give 27 eligible method–condition rows. Each primary macro average weights the eight holo systems equally and retains execution failures as zero workflow successes. A failed program supplies no evidence about biological binding or whether a site exists.

Success means a candidate center lies within 4 Å of any heavy atom of a modeled reference ligand. Candidate ranks are merged at 4 Å using receptor-derived capsid symmetry, with budgets of one and five distinct classes. All native symmetry copies represent one capsid system. The surface baseline is the mean of 100 area-weighted draws with five separated patches per draw; these are conditional Monte Carlo null replicates, not biological observations.

| Method | Top1 success | Top5 success | Top1 common interface recall | Top5 common union recall |
|---|---:|---:|---:|---:|
| fpocket | 37.50% | 62.50% | 41.84% | 61.98% |
| P2Rank | 37.50% | 50.00% | 34.67% | 41.96% |
| Surface baseline | 4.12% | 20.75% | 5.90% | 29.72% |

Common interface recall uses protein heavy atoms within 8 Å of the candidate center and native protein–ligand contacts strictly within 4 Å. Top1 distance and coverage refer to the same closest native copy. For Top5, the union of symmetry-expanded candidate footprints is compared separately against each single native copy, and the maximum recall is reported. Method-specific footprints remain in the detailed data; their differing definitions are not pooled into the common comparison.

## Per-condition outcomes

Each success cell gives Top1 / Top5 at 4 Å. Baseline cells are replicate means.

| Condition | fpocket | P2Rank | Surface baseline |
|---|---:|---:|---:|
| aav9_holo | failed (0 / 0) | 0.00% / 0.00% | 1.00% / 8.00% |
| aav9_apo | failed (0 / 0) | 0.00% / 0.00% | 2.00% / 6.00% |
| 154L | 100.00% / 100.00% | 0.00% / 100.00% | 5.00% / 25.00% |
| 1GNY | 100.00% / 100.00% | 100.00% / 100.00% | 4.00% / 28.00% |
| 1KJL | 100.00% / 100.00% | 100.00% / 100.00% | 5.00% / 19.00% |
| 1QFO | 0.00% / 0.00% | 0.00% / 0.00% | 6.00% / 16.00% |
| 2J1T | 0.00% / 100.00% | 0.00% / 0.00% | 2.00% / 20.00% |
| 1PWB | 0.00% / 0.00% | 0.00% / 0.00% | 0.00% / 16.00% |
| 2VUZ | 0.00% / 100.00% | 100.00% / 100.00% | 10.00% / 34.00% |

P2Rank returned zero candidates for 1QFO and 2J1T; those completed executions remain zero-success observations. Reference quality and assembly limitations apply to every row. No lower-quality case was removed after seeing recovery.

## Distance and candidate-merge sensitivity

The following macro averages retain all eight holo systems. Each cell gives Top1 / Top5 success.

| Method | 3 Å center threshold | 4 Å center threshold | 5 Å center threshold |
|---|---:|---:|---:|
| fpocket | 37.50% / 50.00% | 37.50% / 62.50% | 37.50% / 62.50% |
| P2Rank | 37.50% / 50.00% | 37.50% / 50.00% | 37.50% / 50.00% |
| Surface baseline | 3.25% / 14.75% | 4.12% / 20.75% | 5.12% / 26.75% |

The frozen 2/6 Å merge sensitivities use the same 4 Å center endpoint and retain the original candidate order.

| Method | Merge 2 Å | Primary merge 4 Å | Merge 6 Å |
|---|---:|---:|---:|
| fpocket | 37.50% / 62.50% | 37.50% / 62.50% | 37.50% / 62.50% |
| P2Rank | 37.50% / 50.00% | 37.50% / 50.00% | 37.50% / 50.00% |
| Surface baseline | 4.12% / 20.75% | 4.12% / 20.75% | 4.12% / 20.75% |

## Paired AAV9 conditions

P2Rank's first matching distinct class at 4 Å is rank 11 for holo and 6 for apo, beyond the primary Top5 budget. The raw counts are 774 and 619, respectively, reduced to 25 and 12 symmetry-aware classes. Apo predictions are aligned to the holo evaluation frame with all matching receptor CA atoms; ligands are absent from that fit. The two states supply a paired robustness observation, without an additional independent system.

## Execution history and development amendments

All 9 initial fpocket attempts failed during mmCIF reading (exit 1). The compatibility retry appended a terminal `#` delimiter required by its bundled reader. Atom identities, coordinates, method parameters and the original attempts were preserved. The corrected runs completed for 7 small receptors; both full capsid runs terminated with return code −5 (SIGTRAP), yielding zero predictions. These are workflow failures. The deposited ordered capsid was never cropped to force completion.

A source-only symmetry audit increased the full-capsid maximum CA displacement guard from 0.05 to 0.1 Å for both states before site endpoint scoring. The apo operators reached 0.073316 Å; the chain-fit guard remained 0.05 Å. The primary 4 Å site endpoint, candidate merge distance and 8 Å footprint remained fixed. Original failures and the separately hashed amendment are retained. This development history does not support a pristine prospective confirmatory claim.

fpocket uses unmodified upstream wall-clock seeding; one corrected scientific execution per condition is retained. P2Rank uses fixed seed 42 and one thread. The methods use released configurations rather than equal CPU budgets. Archived raw outputs support result reconstruction without claiming bitwise equivalence of a new native execution.

## Scope of inference

Every predictor received the full supplied protein-only heavy-atom receptor, with all glycans, waters and metals removed uniformly. AAV9 used the full 60-subunit deposited ordered capsid. The seven external receptors are full released files, which can be isolated domains or incomplete biological assemblies. Holo structures retain binding-induced conformational information. Metal removal particularly limits interpretation for calcium-dependent lectins.

The [reference audit](SITE_REFERENCE_AUDIT.md) records weak local density, incomplete glycans, chemistry discrepancies and assembly caveats. Coverage concerns modeled residues only. The development collection and P2Rank training overlap have not passed independent sequence/fold and calibration audits. No family-held-out estimate, experimental nonbinding annotation for other patches, atomic-pose accuracy or ligand specificity follows from these site scores.

## Reproduction

`results/site_discovery/evaluation.json` contains all 27 rows and hashes their detailed/reference records. Each condition retains normalized outputs, original and revised fpocket attempts, raw-artifact archives, frozen input archives and baseline draws. The raw-prediction verifier reparses archived native output; the experiment verifier recomputes reference geometry and scores. This report separately checks source digests, denominator identities and every displayed aggregate against per-replicate evidence.

```sh
python scripts/regenerate_site_report.py --check
```

Evaluated artifact SHA256: `d5043c7160598b100b9c8b78d792ec26a3abde60b88cfa1d288d9fb8ebc57a5f`.
