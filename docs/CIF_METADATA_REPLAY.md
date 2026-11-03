# Replay the CIF metadata amendment

The six archived protein CIF inputs lacked explicit polymer entity/sequence
relationships and a terminal atom-loop comment delimiter used by the pinned
native readers. `heldout_serialization.py` supplies that metadata while
preserving the original author atom identities and all physical atom fields.
Its sequence metadata enumerates **observed residues only**; it does not infer
the full construct or insert unresolved residues.

From a checkout with the repository dependencies installed, run:

```sh
python scripts/verify_serialization_archives.py
python -m pytest -q tests/test_heldout_serialization.py tests/test_serialization_archive_verification.py
```

The verifier requires the repository-pinned `gemmi==0.7.5`. It needs no network,
native reader installation, `tools/`, or `data/processed/` directory. To verify
another copy of this bundle, supply its repository root explicitly:

```sh
python scripts/verify_serialization_archives.py --root /path/to/copy
```

The two immutable preparation archives and their external receipts live under
`results/structural_test/` and `results/structural_test_v2/`. The verifier pins
the nine published source/test/preparation files by SHA256, then checks both
archive hashes and every member's exact name, size and hash. Only ordinary
regular members with canonical relative paths are accepted; duplicate,
missing, unexpected, linked and traversal members fail. Archive members are
read in memory without extraction onto the filesystem.

| Checked item | Count |
|---|---:|
| Original archive members | 99 |
| Amended archive members | 107 |
| Original members with metadata changes | 25 |
| Original members retained byte for byte | 74 |
| Unchanged PDB, PDBQT and SDF members | 34 |
| Unchanged chemistry and canonical ligand-map members | 12 |
| Protein CIFs reproduced exactly by fresh finalization | 6 |

The amendment receipt must enumerate all 99 original members with their
correct old/new hashes. Per-case artifact inventories and aggregate receipts
must agree with the archive contents. The original and amended chemical,
preparation-support and reference-scope decisions remain equal. The extra
eight members are six serialization receipts, the amendment receipt and its
historical configuration proposal.

For each protein CIF, verification compares every original and amended
`atom_site` value by complete author atom identity. Only `label_asym_id`,
`label_entity_id` and `label_seq_id` may change. Elements, coordinates,
occupancies, B factors and every other atom column must be exactly equal;
there is no coordinate tolerance. Selenium-containing MSE retains its literal
residue and element identities.

The verifier also runs the pinned finalizer on each original CIF in a fresh
temporary directory and requires its output bytes to equal the archived V2
CIF exactly. The new proof must equal the historical proof except for its two
input/output path strings. Historical paths are retained as provenance and
need not exist locally. No source receipt is rewritten, and the verifier does
not call `build_source_v2`.

Successful output is JSON containing archive, member, case and atom counts,
including counts of the permitted label-field changes. The checks concern
serialization and immutable source preservation. They do not run predictors,
docking or pose/site accuracy calculations. Archived ligand controls are
preserved and checked for byte identity; no new control result is produced.
