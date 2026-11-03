# Complete development execution archives

These two archives retain the complete ordinary-file trees from the original
156 scheduled runs and the nine corrected neighborhood controls. They include
preparation inputs, raw outputs, failed and partial attempts, recovery shards,
interrupted-attempt evidence, frozen source snapshots and available resource
records. Missing historical resource measurements remain unknown.

The original workflow finished with 153 successful runs and three failures at
its frozen one-hour timeout. The four interrupted earlier attempts remain in
the archive. The corrected controls contain three boxes and three seeds per
box. These overlapping development controls add no independent biological
replicates. AAV remains site-only and has no qualified atomic reference here.

Each manifest covers every regular tar member with its repository-relative
path, size and SHA256, and authenticates the compressed archive. Hydration uses
an empty target prefix. Identical overlap with already published input files
must be verified; conflicting bytes must fail. The exact future derived
`results/blind_pose/qualified_pose_summary.json` is excluded from the raw tree
inventory so generated results remain a separate product.

`package_outputs.py` checks the completed workflow and its raw-evidence audit,
then rejects source-tree changes during packaging. It calls no empirical
reader or native executable. `deterministic_rebuild_verification.json` records
a separate local byte-identical rebuild of both archives, their manifests and
the packaging receipt. This is post-execution packaging, not a prospective
scientific freeze or contemporaneous producer checksum.

The published development workspace plan supplies the exact source commit and
archive pins:

```sh
python scripts/restore_replay_workspace.py --repository . --plan docs/examples/development_workspace_plan.json --destination /tmp/aav-qualified-development-replay
python /tmp/aav-qualified-development-replay/scripts/regenerate_blind_pose_report.py
python scripts/compare_development_report_replay.py results/blind_pose/qualified_pose_summary.json /tmp/aav-qualified-development-replay/results/blind_pose/qualified_pose_summary.json
```

The report needs the project's pinned Python dependencies. It requires no
installed docking program or processed input tree. The separate new-six-case
study retains its own four-workflow completion gate.
