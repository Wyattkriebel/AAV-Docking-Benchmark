# Structural raw-output packaging proposal

This source-only proposal preserves final frozen structural execution evidence. It does not certify terminal execution, inventory or archive real outputs, run scientific software, or authorize accuracy. The new structural accuracy release remains subject to its independent global gate. Root must review [the packager](../results/structural_pose_output_foundation/package_outputs.py), synthetic tests, exact [plan proposal](../results/structural_pose_output_foundation/packaging_plan.proposal.json) and a future terminal-validation receipt before approving an archive run. No readiness or approval receipt has been created.

| Disjoint source root | Frozen scope |
|---|---|
| `results/structural_pose_v3` | Main 216 slots plus nested common-pose fixed-score pool |
| `results/structural_assembly_pose_v3` | 36 additional assembly-control seed slots |
| `results/haddock_paired_cohort` | 108 identities, including failures and first-job recovery history |

The union is nested under main; do not archive it as an overlapping fourth root. The assembly report's 72 comparisons include 36 existing main-comparator slots, rather than 72 new executions.

## Complete physical evidence

Include every ordinary file under these roots except the four named derived-publication exclusions below. Preserve configurations, preparations, archived inputs, commands, manifests, physical raw models, partials, stdout/stderr, AIRs, PSFs, parameters, protocols, resources, wrapper failures, source/control/recovery history, fixed-score inputs and per-pose outputs, progress snapshots and original first-job archive. Preserve failed attempts. Do not reconstruct deleted temporary files or invent missing return codes, CPU time or external wall time.

Installed software, source trees and binaries are outside this scope. External source-contract pins bind prerequisites rather than copying an installation. The published source/protocol foundation is commit `f68775d62551fb76dbd45c98210520b8c752c195`. The proposal binds 41 source/protocol files, including preparation manifests, HADDOCK preflight artifacts, normalized first-job inventory, configurations and plans. It also pins replay-comparison freeze `d31a9018e3b2514d41da4e71959941b4bf7767b14547ff80721b7f173ae29941`. Recheck exact pins before and after packaging; no frozen source or historical manifest changes.

This raw-output increment complements the verified [source/protocol foundation](METHOD_PROTOCOL_SOURCE_REPLAY.md). Exact condition inputs already reside in its 698-member preparation archive and HADDOCK's 215-member preparation and 550-member preflight archives. Restore and verify that foundation first; do not recover inputs from ignored live `data/processed` or installed `tools` trees. The raw archive preserves all currently physical output-root files and the exact historical 240-member overlap, while the published input/source layer supplies prerequisites outside those roots.

## Metadata-only terminal prerequisites

A separate root-reviewed validator must certify readiness without scientific report readers, coordinate/energy interpretation, accuracy, native launches or source writes. File existence alone cannot establish completeness.

- Account for main 216 slots, including actual modern condition-level `execution.json` records and source failures, plus all 108 Vina-Carb records. Missing evidence stays pending.
- Account for 12 assembly conditions with three seeds each, including native and preparation failures.
- Account for all 108 HADDOCK identities through schema-2 progress. The 39 ready native jobs require actual `execution.json` or `wrapper_failure.json`. The other 69 are progress-only preparation/coverage failures; do not fabricate directories or require nonexistent per-job receipts. Unknown original recovery return-code and CPU fields stay null.
- Account for all 216 union source slots and 108 fixed-score pools. The final retained pose-ID set must equal fixed-score record keys, with every record terminal. Count equality alone is insufficient. Historical 1,279-pose progress is not the final denominator.
- Stop every writer to these roots. Require no pending work and absent `results/haddock_paired_cohort/orchestration_recovery/coordinator.lock`. Snapshot stability does not itself prove a stopped writer.

The future receipt has schema 1, status `all_structural_sources_terminal_metadata_only`, `writers_stopped: true`, `pending: []`, `accuracy_computed: false`, `numeric_pose_or_energy_payloads_interpreted: false`, and the full included-plus-excluded snapshot hash. Its exact integer counts are main 216, assembly 36, HADDOCK 108, union_sources 216 and pools 108; retained and recorded are equal nonnegative integers. Root's validator must additionally establish identity-set equality. Booleans and floats are rejected as counts.

## Named exclusions

Only these paths are excluded. Every present excluded file still appears with path, size, hash and reason in the opaque snapshot.

| Path under `results/structural_pose_v3/` | Separate publication |
|---|---|
| `global_reports/` | Future byte-exact four-report packet |
| `global_report_release/generation_receipt.json` | Original generation provenance with that packet |
| `global_report_release/released_reports_binding.json` | Original byte binding with that packet |
| `global_report_release/atom_order_audit/` | Independently reviewed QA publication allowlist |

No cache exclusion is declared. A metadata-only traversal on 2026-11-03 observed no symlinks or special files; actual packaging must recheck. New exclusions require explicit plan review, never implicit globs. Protected supplemental files and unused worksharing proposals are outside these roots and are not traversed.

## Deterministic archives

Sort complete repository-relative POSIX paths. Greedily start a chunk before it would exceed 70,000,000 raw bytes or 1,000 files; reject any individually oversized file. Each gzip must be smaller than 100,000,000 bytes. Raw limits exclude tar headers/padding, so compressed limits are checked separately.

Tar members are regular files with normalized owner, mode and timestamps; gzip has fixed time and no filename. Reject symlink roots/ancestors/descendants, special files, unsafe names, case/path collisions and installed binary content. Equal contents at different paths remain separate evidence. Each complete manifest gives exact tar-relative path, size and SHA256. Members are already repository-relative, so hydration uses an empty target prefix.

Take an included-plus-excluded opaque snapshot before packaging. Verify bytes against it during reads, then require exact complete inventory equality afterward. Recheck source-contract pins and absent locks. Added, removed, renamed or changed included or excluded files abort. The utility never interprets pose, coordinate or score payloads.

Use only the new `results/structural_pose_output_foundation/raw` destination, outside all source roots. Stage privately, create files exclusively and publish the receipt last. An existing destination fails. Cleanup removes only this invocation's unchanged files and preserves foreign replacements. The final receipt binds approval, plan, utility, external terminal validation, snapshots, exclusions, and every archive/manifest hash. Independent reread and private hydration remain publication prerequisites.

## Historical first-job preservation

Keep original archive `c15799f6e92d38cf98d4d4c509a367ce14633faeb6a434dd5d0156ba23895ed7` (3,254,709 bytes), manifest `7bbae28624ea8556e6b903ef638598e996f7fdb6e4695d7e9e4aff49af043ac2` (55,318 bytes), and all 240 members exact. Archive and manifest remain files within the full HADDOCK root. Every historical member must match its current physical counterpart's hash/size. Original manifest rows already have full repository-relative paths; only old tar-relative names need the `results/haddock_paired_cohort` prefix. Never double-prefix or regenerate the original fixture.

## Approval and replay

A future root approval has schema 1, status `approved_terminal_structural_raw_packaging`, `packaging_authorized: true`, exact plan/terminal-validation pins and the reviewed packager hash. It must follow review and bind the actual stopped-writer snapshot. The utility uses the standard library and requires the repository's Python 3.11 or later. Actual invocation after approval is `python3.11 results/structural_pose_output_foundation/package_outputs.py --root REPOSITORY --plan PLAN --approval APPROVAL`, with concrete reviewed files replacing placeholders. The CLI has no synthetic-mode switch.

Before a raw-output PR, independently verify compressed/tar integrity, normalized metadata, exact path/hash/size coverage, caps, snapshots, all 240 historical overlaps, exclusions and fresh private hydration. Synthetic tests use private TEMP fixtures only; run `python3.11 -m pytest -q tests/test_structural_output_packager.py`. No utility fetches, installs or calls scientific readers.

Later replay has two layers: recompute four reports from raw evidence with only explicitly approved derived-roundoff fields; then replay contrast statistics on the byte-exact original four-report packet and generation receipt/binding. No hash tolerance, frozen-source edit or accuracy authorization is introduced here.

Proposal SHA256: `961e5537b1a48223a35cad10f61720db2c0f08772326a1b6f8892d90874bc600`.
