# Independent implementation and evidence review

Reviewed 29 October 2026 (America/Los_Angeles), separately from primary run
orchestration and report authoring. The reviewer also contributed metric controls
and shard validation; this is an internal cross-check, not an external audit.
This review covers reproducibility and the
consistency of bounded development evidence. **The original scientific purpose
remains incomplete.** Passing software checks does not complete its structural,
blind-site, assay-specificity or confirmatory evidence requirements.

## Checks performed

- At the initial final-review checkpoint, **119 tests passed**. Controls cover
  strict atom maps, receptor-frame RMSD, proper rigid transformations, retained
  pose ranking, fixed holo contact denominators, failed runs, stale inputs,
  incomplete manifests, preparation archives, native scoring and report checks.
- Restored the actual AAV9 prepared-input archive into a temporary directory,
  copied all 12 raw docking outputs and frozen manifests/configuration, and
  reran analysis. Both `summary.json` and `pose_metrics.csv` reproduced
  **byte-for-byte**. The committed results were not overwritten.
- Independently hashed all 12 pose files under `results/clean_environment/aav9/`.
  Every file matches the corresponding primary artifact and the hashes in
  `clean_environment_record.json`; its source-manifest hash also matches the
  current frozen source manifest. The record describes a fresh isolated
  environment, source restoration, preparation and 12-run replay on the same
  host. This is implementation reproduction, not independent biological data
  or an external laboratory validation.
- Checked the separate one-seed replay and explicit independent arithmetic in
  `scripts/reproduce_aav9_one_seed.py`. All 20 retained poses reproduce exactly;
  native and translation controls, 35 native contacts, the deposited 12-atom
  GAL map and the 31,080-Cα apo-to-holo transform checks pass. See
  [the reproduction report](INDEPENDENT_REPRODUCTION.md).
- Compared every retained AAV9 pose against all **60 complete deposited GAL
  instances** using the same exact atom-name map and receptor frame. For each
  of the 12 runs, the best RMSD is unchanged from the selected GAL A801
  reference. A different deposited capsid copy therefore does not change these
  particular local-search recovery results. This does not replace a general
  symmetry-aware confirmatory evaluation.
- Independently restored and recomputed **all 15 final structural/assay tables**
  through `scripts.verify_archived_results.py:verify`: every table matched its
  archived counterpart exactly. `python scripts/regenerate_report.py --check`
  also passed for the final generated report. It contains **48/48 completed
  docking runs and 958 retained poses**; four native-coordinate score controls
  remain separate, and the failed 1QFO preparation remains visible.

## Cross-platform verification checkpoint

The subsequent Linux CI run identified a floating-point representation difference
in the 154L per-pose CSV; its JSON summary still reproduced exactly. The default
verifier continues to require exact file bytes. Linux CI now explicitly permits
an absolute tolerance of `1e-10` only for a named whitelist of calculated floating
metrics. CSV identity/count columns and copied values remain exact strings;
JSON booleans, integers, types, shapes, copied assay measurements, resource
records and raw scores remain exact. Input and pose checksums are unchanged.

Review caught and corrected an initial generic numeric coercion that could have
treated numeric-looking identities or copied floats as tolerant metrics. New
controls cover those cases. Independent local validation after that correction
passed **127 tests**, and all **15 tables remain byte-exact locally**, with a
maximum numeric difference of zero under the default verifier.

Linux CI for commit `d530f10` (run `40999522183`) also passed: 127 tests, all
15 tables within the explicit policy, and the generated report current. Its
reported maximum absolute derived-field difference was **4.44 × 10⁻¹⁶**. The
verification record distinguishes byte-exact files from tolerant comparisons
and records actual differences. This numerical portability correction changes
neither scientific inputs/search settings nor outcome decisions or claims.

## Completed 2VUZ execution review

The final 2VUZ manifest contains exactly the three configured seeds, each with
20 complete 90-heavy-atom poses. Pose, log and archived source-manifest hashes
match. Actual commands preserve the frozen two-CPU, exhaustiveness-16,
20-mode, 0.5 Å separation, 10 kcal/mol output-window settings and all original
box coordinates. Independent parsing and explicit receptor-frame RMSD arithmetic
match the final summary for all three seeds:

| Seed | Top1 RMSD Å | Best retained RMSD Å | Best retained rank |
| --- | ---: | ---: | ---: |
| 104729 | 7.2796013817912995 | 5.336553711067763 | 3 |
| 130363 | 7.6444112446791985 | 5.2726146096810655 | 5 |
| 155921 | 6.467147751005333 | 5.656201452093681 | 12 |

The original first seed was preserved. Two separate workers completed the
remaining frozen seeds; this changed scheduling only. The cancelled duplicate
serial attempt remains in `previous_attempts`, its partial log is preserved,
and its unavailable completion timings are explicitly missing. Resource totals
describe final attempts; they do not invent the cancelled work's cost or treat
summed overlapping wall times as elapsed study duration. The seven-residue
resolved glycan remains a fragment reference, not complete nonasaccharide truth.

## Scientific interpretation

The report consistently labels all executed searches as true-site local
development diagnostics. It retains the seven-case external preparation
denominator, including the 1QFO reference-stereochemistry failure. It distinguishes
complete AAV neighborhoods from deliberately incomplete monomer controls and
labels external clash evaluation as limited to the supplied receptor. Source
reuse does not certify external biological assemblies or local density.

The native-coordinate score-only control preserves experimental heavy atoms,
records modeled hydrogen placement and compares the same Vina receptor/box
conditions. It is explicitly excluded from docking success and execution totals.
The native-versus-returned score differences support a prepared-score diagnostic;
they do not identify internal sampling failure or measure binding affinity.

The Gal3 analysis uses three measured binders within one ITC context, reports
alpha/beta reducing-end sensitivity separately, discloses full-length assay
versus isolated-CRD docking, and shows that heavy-atom count already perfectly
ranks this small panel. It does not invent nonbinder labels, convert Vina energy
to KD, or claim discrimination beyond size.

An independent primary-source check confirmed the source-specific claims:

- [Sindrewicz et al. 2019](https://www.nature.com/articles/s41598-019-47658-8)
  supports the ledger's LacNAc 24.0 ± 3.6 µM, lactose 93.5 ± 9.1 µM and
  galactose 3.1 ± 1.1 mM ITC estimates and the disclosed assay construct.
- The [2024 HADDOCK paper](https://doi.org/10.1021/acs.jcim.4c01372) and
  [pinned source repository](https://github.com/haddocking/protein-glycans/tree/344e4a41a300b68b7ee3585895a460c863545ddb)
  support the released-example provenance. Reusing its inputs with Vina is
  correctly distinguished from reproducing published HADDOCK performance.
- [ProtGlycanDock v1](https://arxiv.org/html/2610.26009v1) supports the stated
  ten-method overlap, structural inventory, chain-count exclusion, holdouts
  and predicted receptor/pocket conditions. The repository does not claim that
  these differences alone prove novelty.

No material inconsistency was found in these bounded scientific interpretations.
Unreviewed reference quality, family independence, trained-method leakage,
prospectively frozen held-out execution, tested alternatives and family-aware
precision remain genuine scientific requirements, as the closure ledger states.

## Release checks identified during review

The reviewer reported two reproducibility issues before release: replay examples
must use fresh output directories or verified resume settings now that primary
results are committed, and scientific `.log` artifacts must be included despite
the generic Python ignore rule. Both are resolved: README examples use fresh
replay destinations, and the result-log exception includes scientific logs.
The reviewer also found that panel docking discarded failed statuses; the CLI
now returns a failure status if any configured run fails, consistent with the
single-case command. Downstream panel analysis rejects failed or incomplete runs.
Report replacement is atomic, with a regression control preserving an existing
valid document if writing its replacement fails.

Final 2VUZ aggregation, archive recomputation and report currency checks have
completed. No remaining material scientific-interpretation or reproducibility
defect was found within this bounded reviewed evidence. These checks do not mark
the original scientific purpose complete; the outstanding scientific requirements
remain explicitly recorded in the closure ledger.
