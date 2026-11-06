# Structural reader atom-order QA

This packet publishes a completed source, chemical-label and synthetic-coordinate audit of the new six-case structural readers. It contains no observed docking recovery, native ligand coordinates, empirical pose metrics or method comparison. The old development reader's separate atom-order defect does not occur in the audited new main, HADDOCK, assembly or exact-union paths. The global completion and accuracy-release gate remains mandatory.

The main and assembly parsers explicitly order pose atoms by the same source-map names used to construct the canonical scored mask and `ParentGeometry`. Vina-Carb restores serial-based aliases and then explicitly permutes to that order. HADDOCK emits ligand coordinates in the canonical-key order of the same source parent map. The union reuses exact main per-pose identities and metrics, changing only pose ranking. The detailed source proof is preserved in [original_audit_proof.md](../results/structural_pose_v3/global_report_release/atom_order_audit/original_audit_proof.md).

The packet preserves the original 57,484-byte QA JSON with SHA-256 `c3284c9854f68bbf697f6195a10fc1da2c51dd75fbbb07b170042534be105f64`, its original script and the independent-review proof byte-for-byte. It checked 12 source-profile/case maps, all 19 prepared main name contracts, and retained 17 nonprepared conditions. Independent review checked 13 ready HADDOCK conditions with 332 aliases and nine prepared assembly references plus nine full comparators. These are label-contract counts, not successful docking or biological observations.

The original script records its historical absolute workspace and `/tmp` output paths. The separately published `audit_atom_order.py` adds `--root` and `--output`, with output required outside the audited root. Choose a fresh external filename; the script preserves the original writer behavior and can replace an existing external path. Only import/CLI/root/output path handling changed; the exact patch and packaging amendment record both script hashes. Every original algorithm statement and assertion remains unchanged. Run without Python's assertion-disabling `-O` option, using the repository's pinned scientific dependencies:

```sh
python results/structural_pose_v3/global_report_release/atom_order_audit/audit_atom_order.py \
  --root . --output /tmp/aav-atom-order-replay.json
cmp results/structural_pose_v3/global_report_release/atom_order_audit/original_audit.json \
  /tmp/aav-atom-order-replay.json
```

The runner reads archived chemical maps and SDF graph symbol/index fields, plus prepared names/maps. It hashes source archives opaquely and ignores source SDF coordinate columns. Runtime parser and graph-order controls use synthetic coordinates only. It does not open actual pose outputs or native ligand reference ledgers, call a scorer/report/readiness probe, fit a protein or compute RMSD, contacts or geometry diagnostics.

`portable_replay_verification.json` records byte-exact replay from a private source/name-map fixture with tools and processed inputs absent and an unrelated working directory. An inside-root output attempt was rejected before any audited input was read. These checks reproduce software QA, not structural performance. `packaging_amendment.json` and the independent publication review distinguish original evidence from the later path-only publication work; none replaces an original scientific freeze, the historical735-file inventory or frozen metric implementation.
