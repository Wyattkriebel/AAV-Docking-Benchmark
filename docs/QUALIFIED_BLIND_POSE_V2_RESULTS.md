# Qualified predicted-site and oracle-site pose results

V2 corrects atom ordering by a complete named coordinate-row permutation. This implementation amendment follows invalid v1 results and corrected diagnostic outcomes; it is not pre-outcome.

This development comparison separates complete modeled ligands from independently qualified modeled cores. No core result is called whole-ligand recovery.

The fixed workflow contains 27 conditions and 156 scheduled seed runs. Atomic accuracy is eligible for six external complexes: two complete models and four modeled cores. AAV9 and 154L remain atomic-reference ineligible. Unsupported 1QFO chemistry remains zero workflow recovery within the eligible core denominator.

Each cell below is the mean of three seed-specific success indicators per complex, then the equal-weight mean across complexes. The primary cutoff is receptor-aligned RMSD ≤2 Å. Site budget means searched candidate boxes; pose rank means energy-ranked returned poses within that site budget.

## Complete modeled ligands

| Method | Site budget | Global Top1 | Global Top5 | Global best returned | Local Top1 | Local Top5 | Local best returned |
|---|---:|---:|---:|---:|---:|---:|---:|
| fpocket | 1 | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% |
| fpocket | 5 | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% |
| p2rank | 1 | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% |
| p2rank | 5 | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% |
| oracle | 1 (oracle) | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% |

## Qualified modeled cores

| Method | Site budget | Global Top1 | Global Top5 | Global best returned | Local Top1 | Local Top5 | Local best returned |
|---|---:|---:|---:|---:|---:|---:|---:|
| fpocket | 1 | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% |
| fpocket | 5 | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% |
| p2rank | 1 | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% |
| p2rank | 5 | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% |
| oracle | 1 (oracle) | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% |

## Per-complex results

The following uses the five-site predicted budget and the one-box oracle control. Each percentage averages all three frozen seeds, retaining failed eligible runs as zero success.

| Complex | Qualified scope | Method | Global Top1 / Top5 / best returned | Local Top1 / Top5 / best returned |
|---|---|---|---|---|
| 1GNY | full_ligand (46 atoms) | fpocket | 0.00% / 0.00% / 0.00% | 0.00% / 0.00% / 0.00% |
| 1GNY | full_ligand (46 atoms) | p2rank | 0.00% / 0.00% / 0.00% | 0.00% / 0.00% / 0.00% |
| 1GNY | full_ligand (46 atoms) | oracle | 0.00% / 0.00% / 0.00% | 0.00% / 0.00% / 0.00% |
| 1KJL | modeled_core (11 atoms) | fpocket | 0.00% / 0.00% / 0.00% | 0.00% / 0.00% / 0.00% |
| 1KJL | modeled_core (11 atoms) | p2rank | 0.00% / 0.00% / 0.00% | 0.00% / 0.00% / 0.00% |
| 1KJL | modeled_core (11 atoms) | oracle | 0.00% / 0.00% / 0.00% | 0.00% / 0.00% / 0.00% |
| 1QFO | modeled_core (31 atoms) | fpocket | 0.00% / 0.00% / 0.00% | 0.00% / 0.00% / 0.00% |
| 1QFO | modeled_core (31 atoms) | p2rank | 0.00% / 0.00% / 0.00% | 0.00% / 0.00% / 0.00% |
| 1QFO | modeled_core (31 atoms) | oracle | 0.00% / 0.00% / 0.00% | 0.00% / 0.00% / 0.00% |
| 2J1T | modeled_core (21 atoms) | fpocket | 0.00% / 0.00% / 0.00% | 0.00% / 0.00% / 0.00% |
| 2J1T | modeled_core (21 atoms) | p2rank | 0.00% / 0.00% / 0.00% | 0.00% / 0.00% / 0.00% |
| 2J1T | modeled_core (21 atoms) | oracle | 0.00% / 0.00% / 0.00% | 0.00% / 0.00% / 0.00% |
| 1PWB | full_ligand (23 atoms) | fpocket | 0.00% / 0.00% / 0.00% | 0.00% / 0.00% / 0.00% |
| 1PWB | full_ligand (23 atoms) | p2rank | 0.00% / 0.00% / 0.00% | 0.00% / 0.00% / 0.00% |
| 1PWB | full_ligand (23 atoms) | oracle | 0.00% / 0.00% / 0.00% | 0.00% / 0.00% / 0.00% |
| 2VUZ | modeled_core (47 atoms) | fpocket | 0.00% / 0.00% / 0.00% | 0.00% / 0.00% / 0.00% |
| 2VUZ | modeled_core (47 atoms) | p2rank | 0.00% / 0.00% / 0.00% | 0.00% / 0.00% / 0.00% |
| 2VUZ | modeled_core (47 atoms) | oracle | 0.00% / 0.00% / 0.00% | 0.00% / 0.00% / 0.00% |

## Interpretation and verification

No independent ligand fit is used. Global alignment uses all corresponding protein backbone N/CA/C/O atoms. Local alignment uses backbone atoms of native qualified-core contact residues, solely during evaluation. Chemically valid full-parent automorphisms preserve source graph/stereochemistry; the same chosen map supplies RMSD and contacts. Native controls give zero RMSD, full native contacts and 5 Å error after the known [3,4,0] Å translation.

Per-pose 3D stereochemistry, bond-length changes, ring-distance changes, linkage distances and full-parent/protein severe overlaps are recorded separately. These diagnostics do not certify ring pucker, glycosidic torsion energies or overall physical validity. Primary RMSD success is not silently converted into a pose-validity claim.

The geometry fields named change_from_input and glycosidic input_distance_A use the archived independent legacy SDF conformer. It shares the executed ligand's initializer lineage but precedes PDBQT and box-translation coordinate rounding; it is not an exact executed-box PDBQT baseline.

All failures retain their eligible workflow denominator. Three seeds, capsid copies and repeated boxes are not independent biological observations. Full-model and core strata are small development samples; no family-level reliability interval or universal docking conclusion is earned.

Inputs restore from checked-in preparation archives, while native coordinates enter only this evaluator after execution completeness is verified. Frozen configuration, reference ledger, original structural sources, input/pose hashes, commands and raw Vina log tables are checked. The summary retains every pooled candidate identity and source checksum.

Machine-readable results: [qualified_pose_v2_summary.json](../results/development_pose_report_release/qualified_pose_v2_summary.json). Atomic eligibility: [atomic_reference_ledger.json](../benchmarks/atomic_reference_ledger.json).
