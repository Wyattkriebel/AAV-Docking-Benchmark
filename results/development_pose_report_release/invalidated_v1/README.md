# Invalidated V1 analysis — audit evidence only

These exact files preserve the first derived development report, its generation
and byte-replay receipts, and the independent checks that rejected it. V1 used
blind-box atom traversal order where its legacy parent graph expected canonical
atom-name order. Its RMSD, contact and geometry values are invalid.

Raw docking files, commands, run counts and failure records remain valid.
The corrected V2 evaluator applies a complete named coordinate-row permutation
under an explicit implementation amendment. Both invalid V1 results and the
independent corrected diagnostics were already known when that amendment was
made; this is a post-outcome implementation correction.

Use the [V2 reproduction guide](../../../docs/QUALIFIED_DEVELOPMENT_REPRODUCTION.md)
for accepted corrected output. The historical files here are retained without
editing their bytes. Their original paths and hashes are recorded in
`../v1_invalidation.json`.
