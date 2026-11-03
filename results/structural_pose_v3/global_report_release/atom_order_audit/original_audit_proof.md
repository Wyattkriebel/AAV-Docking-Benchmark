# New six-case atom-order audit — source and labels only

Scope: the frozen new paired-engine, HADDOCK, assembly-control and exact-union readers, their shared heldout scorer/ParentGeometry and fixed-coordinate export path. No actual pose/native ligand coordinates or reference ledger were read. No workflow probe, report generator, RMSD, protein fit, contact computation, geometry diagnostic or native executable was called. Prepared chemical name/maps and SDF graph symbol/index fields only were read; coordinate columns were ignored. All coordinates used in runtime controls were synthetic. No repository files changed.

## Finding

The frozen new readers explicitly establish a single per-source atom ordering before the common metric call. They do not use the old development reference as a second independent order. No analogous atom-order defect was found in the audited new reader paths. This clears this specific source/input-label concern; it is not a report of docking accuracy, a substitute for the global completion guard or a universal correctness claim.

## Proof through the shared scorer

Let S be the ordered atom list in the profile archive's `<case>/ligand_atom_map.json`, and N be its ordered `name` list.

1. `heldout_pose_report._source_context` at459–463 constructs both `reference.names=N` and `reference.canonical_atom_map=S`; it derives `mapping` and passes that SAME N/S to `ParentGeometry`. There is no call to legacy development reference construction.
2. `canonical_mapping` at114–147 builds `canonical_atom_key -> names.index(atom.name)`. Its mask is the ordered qualified `scored_atom_map` key list. Every native-copy key list must equal that mask in order. Therefore `mapping.indices` indexes N in exactly the native scored-mask order, independently of native/source serialization order.
3. `ParentGeometry` at395–417 maps N to RDKit heavy order through each S record's `rdkit_index`, stores the inverse permutation with `argsort`, and expresses every scored chemical automorphism back in N indices. `diagnostics` at428 uses that inverse permutation before graph geometry. The shared new `score_pool` at192 and215–216 indexes the parsed N-ordered array with those N-indexed maps. Contacts at222 use the same indexed order.
4. The actual12 strict/repair archive maps have unique complete name/key/RDKit-index bijections. All12 match source SDF heavy-index/element graph fields, without reading coordinate columns. Strict/repair maps are identical per case. This metadata check did not invoke `_source_context`, `canonical_mapping` or ParentGeometry with actual study input/reference coordinates.

## Main paired modern Vina and Vina-Carb

- `generate_paired_report` at526–536 obtains S from the pinned profile source archive and requires normalized prepared canonical map equality as a LIST, preserving order, not only set equality. It also binds the initializer SDF to the same archived graph.
- `read_execution_pair` at420–430 takes N in that prepared list order and supplies it to both native-output readers.
- `read_modern` at344–347 requires the prepared `reference.names` to equal N exactly. `read_poses` at400 uses N; `engine.read_poses` at51–53 calls `metrics.mapped_coordinates`, which at35–49 rejects missing/extra/duplicate names and explicitly returns coordinates in required N order regardless of PDBQT record/torsion-tree order.
- `additional_method.read_carb_poses` at163–194 restores generated names through the invariant serial map, returning serial order. Crucially `canonical_carb_poses` at166–178 then validates the complete name universe and explicitly permutes that array into N. The serial order is never assumed to equal N.
- Main contexts at558 use the same profile archive S; the resulting parsed arrays are passed unchanged to new shared `score_pool` at570–573. No legacy reference is substituted.
- Actual source/name check: all19 prepared references and prepared canonical maps equal the same archive S/N order. The other17 of36 source conditions have no prepared pose order and remain explicit failure/unavailable cases.

## General heldout reader

`heldout_pose_report.read_workflows` at393–395 requires each prepared `reference.canonical_atom_map == S` and `reference.names == N`; at442 it parses poses using those reference names. At637–646 the metric context is built from that SAME profile source archive. Thus this direct heldout path also avoids the legacy mismatch even if physical PDBQT atom records are rearranged.

## HADDOCK aliases and source parent order

`haddock_pose_report._condition_map` at630–656 binds source-only preflight aliases to the independent parent graph, requires every canonical heavy key exactly once and preserves full-parent chemistry. `parse_model` at204–220 takes ordered `parent_keys` directly from the passed parent atom list and creates bijections between CNS identity and canonical key. Its return at261 is explicitly `[ligand[by_canonical[k]] for k in parent_keys]`, not raw PDB record order, CNS name order or sorted aliases. In `generate_report`, source contexts at737 and `parent` at743 come from the SAME `sources[profile][case/ligand_atom_map.json]`; at762 that parent list controls parsing, and at767–769 the array is passed to the context from the same S. CNS chain/residue aliases affect lookup only, never the final array order.

Independent archive-alias and synthetic reordered-HADDOCK checks are recorded in the accompanying reviewer handoff; no actual native model was parsed in this audit.

## Assembly control

`assembly_control_report` at241–253 loads pinned monomer and full source profiles and validates paired ligand sources. At303–310 each row passes `[a.name for a in S]` to `read_modern`, which enforces the prepared reference-name equality and performs explicit parser mapping. At325–338 contexts and scored arrays are selected using the SAME profile/case source. The alternate observed full-assembly protein evaluation map changes protein truth only; it does not replace ligand S/N or create a legacy order.

## Exact union metrics and rescoring

- The union index uses `read_execution_pair` (`paired_common_pose_rescoring` at86–100), so its source model enumeration has the same validated modern/VC name contract. The union index stores identities and raw source/MODEL hashes, not a new ligand coordinate order.
- Source model identity is condition/seed/source-algorithm/original-model-rank. `original_source_orders` records those original IDs; sorting the global index by `pose_id` does not change original model-rank identity or pose atoms.
- VC raw alias restoration at166–204 maps by serial and copies only generic identity fields from prepared atoms. It checks exact serial universe, generated names/types/charges and coordinate preservation. `export_source_pose` at212–237 rechecks input/export/raw/model hashes, restores VC aliases if required, and calls the unchanged `export_vina_carb_pose`.
- `additional_method.export_vina_carb_pose` at225–270 reconstructs the fixed source torsion tree by serial and verifies each serial's name/residue/type/charge before transferring its coordinate columns. The fixed scorer does not use a native RMSD mask.
- `paired_common_pose_report.join_main_pose_metrics` at188–248 requires exact condition/source/model identity and main/union pose universes, canonical-map/input hashes, record/raw/model hashes, candidate/rank/native-score identity. It copies the existing fresh main per-pose metric record by exact pose ID; it does not fit or remap coordinates.
- `summarize_ranked_pools` at265–281 ranks those IDs using fixed coefficients and changes only pooled rank and ranking score. It cannot mix atom order because it never rebuilds ligand coordinate arrays. Best-retained invariance is separately asserted at282 onward.

## Runtime controls and reproducible metadata evidence

`/tmp/aav-new6-atom-order-audit.py` writes `/tmp/aav-new6-atom-order-audit.json` (SHA256 `c3284c9854f68bbf697f6195a10fc1da2c51dd75fbbb07b170042534be105f64`). It records all source code pins,12 source name/key/index maps, source graph index checks,19 prepared-name contracts and17 retained nonprepared conditions. It never opens an actual pose file or native reference ledger.

Its synthetic5-heavy-atom graph deliberately uses distinct source-name, RDKit, raw-record and VC-serial orders. The qualified mask is nontrivial `[2,4,0]`. Controls verify:

- modern parser returns source-name order despite scrambled raw records;
- VC parser returns a different serial order, then the explicit canonical wrapper restores source-name order;
- graph name/heavy inverse permutation reconstructs the synthetic RDKit coordinates and identity scored automorphism uses the same mask;
- VC-to-generic union alias roundtrip preserves every synthetic atom coordinate and produces the same N-ordered parsed array;
- no actual RMSD, fit, contact, geometry diagnostic, report or workflow probe is called.

No frozen implementation amendment is indicated for this particular new-six-case atom-order concern. The old development c773 mismatch remains a distinct defect requiring its separately versioned fix and honest invalidation/amendment; these new-path checks do not validate the old output.

## Independent source check

`/root/additional_docking_method/contrast_contract` independently confirmed the HADD/assembly paths and reviewed the synthetic/metadata script. All13 ready HADDOCK conditions/332 alias atoms have bijective CNS identities and canonical parent keys; RDKit indices, elements, formal charges, residue identities and stereo fields match the source parent maps. All9 prepared assembly references plus9 prepared full comparators have exact source name order. Both assembly cases' maps match across strict/repair/full/monomer profiles;16 opaque SDF/PDBQT hashes confirm identical paired ligand bytes. Synthetic controls shuffle CNS aliases, PDB/PDBQT records and parent-name order, and recover full/partial masks through parser-to-ParentGeometry composition. The reviewer found no order defect or need for a versioned new-reader amendment, and called no scoring, fit, RMSD, contacts or geometry diagnostics.

Status: this specific new-six-case atom-order concern is cleared. The original global completeness/accuracy authorization remains mandatory. This clearance does not rehabilitate legacy c773 development metrics or certify any measured study outcome.
