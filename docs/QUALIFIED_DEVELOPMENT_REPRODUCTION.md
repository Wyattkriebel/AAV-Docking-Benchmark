# Corrected development V2 replay

The corrected OLD V2 report uses a bijective atom-name permutation to put raw
blind-box coordinates in canonical parent order before evaluation. The V1
implementation compared incompatible atom orders. Its invalidated results,
release receipt and workspace plan remain under
[`invalidated_v1`](../results/development_pose_report_release/invalidated_v1/)
for audit; they are not the primary development result. The raw archived docking
executions remain valid inputs.

This is a post-outcome implementation correction. Both the invalid V1 outcomes
and corrected diagnostic outcomes had already been seen before the V2 amendment
and freeze. The amendment records that timing; the corrected development report
is not a prospective or confirmatory study.

Once the corrected source/data commit, freeze and canonical outputs are
published, reproduce the report from the repository root with:

```sh
python scripts/verify_qualified_development_replay.py
```

The default requires identical JSON bytes. An explicit cross-platform comparison
allows an absolute difference of at most `1e-10` only at the named derived
metric, fit, geometry, contact-fraction and fractional-statistic paths in the
[frozen comparison policy](DEVELOPMENT_REPORT_REPLAY.md):

```sh
python scripts/verify_qualified_development_replay.py --numeric-tolerance 1e-10
```

Native scores, copied coordinates, atom maps/order provenance, identities,
counts, thresholds, success decisions, source hashes and amendment provenance
remain exact. Markdown must always match byte for byte. The frozen comparator's
schema-source identifier still names V1 because its derived field paths are
unchanged; the V2 evaluator and new protected provenance are bound separately
by the V2 source freeze. The verifier does not widen or mutate that policy.

The verifier checks the literal SHA256 pin for
[`evaluation_v2_freeze.json`](../results/development_pose_report_release/evaluation_v2_freeze.json),
every frozen source file and the checksum-pinned implementation amendment. It
checks the canonical JSON and Markdown against the two output pins in
[`generation_v2_receipt.json`](../results/development_pose_report_release/generation_v2_receipt.json),
and requires a successful atomic generation tied to the exact workspace-plan
hash and source commit.

It then exports that commit into a fresh temporary directory and verifies and
hydrates exactly the OLD156 and corrected-nine output archives declared in
[`development_v2_workspace_plan.json`](examples/development_v2_workspace_plan.json).
Every archive and member is checked against its manifest. No installed `tools`
or `data/processed` tree is copied. The source commit must contain the V2 freeze,
evaluator and amendment, and must not already contain either derived V2 output.
The verifier invokes the exported V2 CLI with both required amendment arguments,
then compares the recomputed files with the canonical
[JSON](../results/development_pose_report_release/qualified_pose_v2_summary.json)
and [Markdown](QUALIFIED_BLIND_POSE_V2_RESULTS.md).

This replays the archived OLD development evaluation; it does not run docking,
curate references or compute the separate new six-case structural study. It
requires the repository's Python dependencies and local Git objects, with no
network downloads or native docking executables. Add `--output /path/to/new.json`
to save the verification receipt; an existing receipt is never overwritten.

The verifier contains a literal pin for the finalized V2 evaluation freeze.
The freeze excludes this verifier from `source_files` to avoid a hash cycle;
the verifier is pinned by its own Git commit. The source/data commit precedes
official V2 generation. The later result release records the workspace plan,
generation receipt and canonical outputs. Full replay is validated only after
those artifacts exist; source-only and synthetic checks do not establish that
the corrected scientific report has reproduced.
