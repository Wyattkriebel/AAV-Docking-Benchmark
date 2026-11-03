# Restore a verified replay workspace

The standalone restore utility exports an explicit local Git commit and
hydrates only checksum-pinned archive members. It uses Python's standard
library and local Git. It does not execute exported programs or call a metric
reader. The live working tree, index and untracked files are read neither as
export inputs nor as archive sources.

Create a destination that does not already exist and is outside the source
repository. The example restores the published preparation snapshot and its
698-member input archive:

```sh
python scripts/restore_replay_workspace.py --repository . --plan docs/examples/preparation_workspace_plan.json --destination /tmp/aav-preparation-replay
python /tmp/aav-preparation-replay/scripts/verify_prepared_inputs.py --root /tmp/aav-preparation-replay
```

The second command needs the project's pinned dependencies. It runs only the
input-correspondence verifier. Neither command runs predictors, docking or
experimental accuracy calculations. The workspace contains no `.git` metadata,
`tools/` or `data/processed/` tree.

Use a complete local repository containing the plan's exact commit and all its
objects. Shallow histories may lack that ancestor. Partial/promisor repositories
are rejected, and the helper never fetches missing objects. Git blob objects
are read directly and checked against their object IDs, so export-ignore,
export-subst and working-tree filters cannot alter the snapshot.

The plan has exactly three fields: `schema_version:1`, a full 40- or 64-digit
`commit`, and an `archives` list. Each entry declares `archive` and `manifest`
pins (`path`, `sha256`, optional `bytes`) plus an explicit `target_prefix`.
Source paths name files inside the Git export. An empty prefix keeps tar member
names as repository-relative destinations; another prefix is prepended to
each member name.

Each coverage manifest must declare `archive_sha256` and a `files` list whose
records contain every raw tar member's `path`, `bytes` and `sha256`. If an
`archive` basename is declared, it must match the pinned source. Links,
directories, sparse members, path traversal, noncanonical names, duplicate or
case-aliasing paths, missing/extra members and size/hash mismatches fail.
Members are streamed to ordinary files without `extract` or `extractall`.
Source Git links and submodules are rejected too. Tools, processed-tree and
reserved `.replay/` paths cannot be exported or hydrated.

Existing files are never replaced. Identical Git files or files supplied by an
earlier verified archive can be checked and reused, with each source retained
in provenance; different bytes or file/directory conflicts fail. Every archive
and manifest must be listed explicitly. The utility does not discover local
untracked archives or infer a scientific completion state.

All validation happens in a private sibling staging directory. Publication
exclusively claims the requested destination and uses exclusive file creation,
with the provenance receipt written last. A destination created by another
process is preserved. Integrity failures leave no destination. If publication
itself fails, cleanup removes only unchanged files created by this call and
preserves unrecognized concurrent changes. A workspace is complete only after
the call succeeds and `.replay/provenance.json` exists.

The receipt records the exact commit/tree, every exported blob ID and file
hash, the plan hash, every archive/manifest pin, every hydrated target hash and
size, and whether each member was created or identically reused. Frozen
readers and scientific plans remain unchanged.

For private test fixtures that already copied explicit source files,
`hydrate_verified_archives(root, entries)` accepts the same archive entries.
The caller must exclusively own a fresh, unversioned temporary source tree.
Git working trees, protected directories and symlinks are rejected. Validation
uses staging and only new files are published; existing bytes remain intact.
Its returned audit labels reuse as `existing_identical` and makes no Git-commit
provenance claim.

Future terminal-output bundles can use the same coverage format after their
completion gates. Older archives that embed an original manifest, or whose
original inventories use another path base, need a separate complete packaging
coverage manifest. List every unchanged member, including the embedded
manifest itself, and supply the appropriate prefix. No frozen original
manifest or numerical evidence needs to be rewritten.
