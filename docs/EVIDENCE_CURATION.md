# Evidence and eligibility ledger

Curated 28 October 2026. The machine-readable ledgers distinguish a downloadable
example from a reviewed experimental reference. An unchecked quality field fails
confirmatory eligibility; it does not silently become a pass.

`benchmarks/reference_ledger.json` records seven nonviral development candidates,
the AAV9 diagnostic, a regional AAV2 evidence case, and a sulfated AAV-DJ candidate.
The seven nonviral cases come from the authors' [released HADDOCK benchmark](https://github.com/haddocking/protein-glycans/tree/344e4a41a300b68b7ee3585895a460c863545ddb),
pinned before new runs. Its prepared bound/unbound files provide a useful external
software reproduction; they do not certify family independence, local density, or
the original apo deposition identity for this study.

`benchmarks/sources/glycan_connectivity.json` records deposited branch bonds,
leaving atoms, author/label residue mappings, source revisions and SHA-256 checksums.
Bonds come from mmCIF declarations, not distances. The external preparation
pipeline must still verify its ordered atom maps and stereochemistry against these
sources. No reference-derived coordinates are eligible as independent input
conformers for a docking claim.

The [2VUZ deposition](https://www.rcsb.org/structure/2VUZ) has seven glycan residues
in its declared graph although its title describes a nonasaccharide. Its modeled
core can support a specifically defined source-reproduction check. Missing tails
cannot be evaluated as exact atomic coordinates. Local density review must define
any scientifically scored core before test docking. [1PWB](https://www.rcsb.org/structure/1PWB)
and 2VUZ also require explicit treatment of recognition-associated calcium;
a dry preparation that deletes calcium does not reproduce the metal-containing
experimental context.


`benchmarks/sources/reference_validation.json` archives machine-extracted current
wwPDB validation diagnostics. [154L validation](https://files.rcsb.org/pub/pdb/validation_reports/54/154l/154l_validation.xml.gz)
contains ligand RSCC values 0.602, 0.794 and 0.668 plus an EDS R mismatch
(0.3099 versus depositor 0.159). [1KJL validation](https://files.rcsb.org/pub/pdb/validation_reports/kj/1kjl/1kjl_validation.xml.gz)
also flags an R mismatch. These warnings require interpretation before a
confirmatory atomic-pose claim. AAV9 and AAV-DJ validation does not provide the
X-ray ligand RSCC fields; missing values are not evidence of good local fit.

The confirmatory 24–36-complex cohort is not frozen. No candidate in this inventory
has been promoted to confirmatory eligibility. Required work includes validation
and local density review, assembly selection, complete chemical mapping,
independent apo provenance, family clustering and leakage review. A deterministic
inventory is useful, but it is not evidence that those checks passed.

## Specificity sources

`benchmarks/assay_ledger.json` includes a three-ligand quantitative affinity panel
for exploratory evaluation and preserves exclusions for unqualified sources. No
experimental nonbinder labels are invented.

[Mietzsch et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC3958061/)
provide an assay lead with synthetic heparin alternatives, including qualitative
serotype-specific binding comparisons. AAV9 was excluded. The numbered compounds
still need complete chemical transcription, numerical measurements and a fixed
label rule. Its heterogeneous natural-heparin control is not a single exact ligand.

Current CFG pages for [AAV2 primscreen_5347](https://www.functionalglycomics.org/glycan-array/1005346)
and [AAV1 primscreen_5777](https://www.functionalglycomics.org/glycan-array/1005776)
both explicitly mark their result nature **Inconclusive**. Their public pages
identify the workbooks, but current downloads require reCAPTCHA. These flags must
be reconciled with any later use of the original measurements. The ledger also
preserves a labeling discrepancy between the record and article description.
Public visibility alone does not qualify binary binding labels.

[Zhang et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC3859860/)
supply SPR/competition context using heterogeneous heparins and size fractions.
These measurements do not automatically define exact-glycan alternatives or a
monovalent-affinity benchmark. None of these sources permits transferring labels
between serotypes or equating docking scores with dissociation constants.

## October 2026 overlap audit

[ProtGlycanDock v1](https://arxiv.org/html/2610.26009v1), posted 19 October 2026,
already evaluates ten methods on 1,781 curated noncovalent complexes with separate
protein and glycan holdouts. Its filtering removes structures with more than four
protein chains. Its Vina arms use predicted receptors and predicted pockets.
Consequently, a general glycan docking leaderboard substantially overlaps this
work. Assembly dependence, controlled information regimes and independently
assayed specificity remain scope distinctions to investigate; they are not proof
of novelty. No universal docking reliability claim or new performance estimate is
supported by this overlap audit.

## A recoverable quantitative panel

[Sindrewicz et al.](https://www.nature.com/articles/s41598-019-47658-8)
report same-assay galectin-3 ITC affinities for LacNAc, lactose and galactose.
`benchmarks/galectin3_itc_panel.json` preserves the measurements and conditions.
All three are binders; this panel supports descriptive affinity ranking, not a
binder/nonbinder classifier. Docking to 1KJL approximates the full-length assayed
protein with its CRD. Free reducing ends require anomer sensitivity, and the source
does not specify exact PBS pH/ionic composition. Three compounds and one target
cannot establish general reliability. Measured rank and heavy-atom count have
the same ordering here, so a successful raw rank alone cannot show information
beyond molecular size.
