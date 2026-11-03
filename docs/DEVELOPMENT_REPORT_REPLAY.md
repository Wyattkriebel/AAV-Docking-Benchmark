# Qualified development report comparison

The old 156-run development evaluation has a separate output schema from the
four newer structural reports. Its read-only comparison command is:

```sh
python scripts/compare_development_report_replay.py saved/qualified_pose_summary.json replay/qualified_pose_summary.json
```

The default validates both JSON files and requires identical bytes. An explicit
`--numeric-tolerance 1e-10` permits absolute roundoff only at reviewed derived
paths; the maximum permitted tolerance is `1e-8`, with no relative term.
`--structural` compares exact parsed types and values while ignoring formatting.
The finite binary64 parsing and exact-ratio arithmetic are documented in
[REPORT_REPLAY.md](REPORT_REPLAY.md).

The separate `qualified_development_derived_paths_v1` policy follows the frozen
`qualified_pose.py` evaluator, SHA256
`c77333c8fbda7fe81f79fbe99211a65f3c91f51382daaa80f2fb89bd12dfc525`.
It covers per-seed global/local pose and rank-summary RMSDs, contact fractions,
bond/ring and linkage distances, receptor-fit diagnostics, fractional seed/case
averages, and descriptive paired contrasts. Raw seed success booleans remain
exact even though their fractional averages are eligible for explicit tolerance.
Native Vina energies, raw atom-mapping coordinates, analytic zero native RMSD,
identities, ranks, maps, counts, thresholds, decisions, measurements and hashes
also remain exact. Unknown paths have no tolerance.

The implementation compares the protected original tree exactly, verifies that
the approved float path sets/types match, then passes only those scalar pairs
to the unchanged strict comparison core. Its internal scalar carrier exists
solely to reuse comparison arithmetic; it computes no scientific quantity.
Diagnostics are mapped back to original development-report paths. Neither the
frozen evaluator nor the newer report comparator's policy is modified.

This command compares two already generated files. It imports only standard
library code and the frozen comparison core, and works with `python -S`. It
does not invoke an empirical evaluator, establish workflow completion or replace
raw-source verification. Generate the development report only after its own
156-run terminal gate has been verified.
