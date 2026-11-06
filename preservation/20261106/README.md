# Local-work preservation — 6 November 2026

This is an archival checkpoint, not completion of the scientific study.
The focused user-facing comparison remains in PR #14. This preservation branch
retains the earlier research implementation and unmodified support proposals.

## Public coverage

- All earlier local commits were already reachable on the fork.
- Ten previously untracked scheduling/opaque-packaging support files are now
  committed unchanged, including their synthetic tests and proposal documents.
- 14,370 existing nonviral per-case files are preserved byte-for-byte in 16
  release assets, totaling 210,377,519 compressed bytes. No scientific method,
  prediction, pose-accuracy calculation or new biological study was run.

The data assets are at the [preservation release](https://github.com/AlecKriebel/AAV-Docking-Benchmark/releases/tag/preservation-2026-11-06).
[receipt.json](receipt.json) records archive sizes and SHA256 values;
[members.ndjson](members.ndjson) records every original path, size and SHA256.
Archives contain repository-relative regular files only, with no symlinks or
path traversal. Download them and extract into a fresh clone to restore these
case artifacts. The release avoids adding the large archives to Git history.

This is an archival selection from declared nonviral case directories. It does
not define a new benchmark cohort, change a frozen scientific estimand, or
qualify the correctness of the historical outputs. Nested archives and files
with mixed-scope markers were kept private. Support proposals retain their
historical warnings and must not be treated as current execution instructions.

## Private coverage and deletion

This is **not a complete public backup**. Of the 16,238 originally untracked
paths, 14,380 are covered by the public data archives or support commit and
1,858 remain private. Ignored files, runtime installations, Git internals,
orphan objects, blocked expansion files and unpublished mixed/viral artifacts
also remain in the separate private backup.

Before publication, both working folders were independently cloned to the
owner's private Repository-Backups directory. All 220,736 nonvolatile regular
files (7,070,109,269 bytes), path sets, modes and symlink targets were checked.
Finder .DS_Store caches were copied but excluded from equality checks because
they change during inspection. The private backup includes the shared .git
directory, including objects omitted by a normal clone or git bundle --all.

Keep the private backup if deleting the original working folders. The focused
working folder shares the main repository's Git directory, so deleting the
main folder alone breaks Git access in that focused checkout. Private restore
instructions and the byte-verification manifest accompany the backup.

No original working folder or research artifact was deleted. The only cleanup
performed was temporary test/publication staging created during preservation.

## Validation

All archive members were reopened and checked against their recorded SHA256.
The ten support files match the captured local bytes. Their 193 synthetic tests
passed with per-test temporary-directory cleanup; the first attempt exhausted
local temporary disk space. Original scientific sources and frozen protocols
are unchanged. Scientific workflows are not rerun by this archival commit.
