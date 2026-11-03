# Atomic-reference qualification of the seven released development cases

On 3 November 2026, a source-only audit qualified **two full modeled ligands and
four explicitly delimited modeled cores**. 154L supplies no atomic reference under
these rules. Every original case remains in the seven-case denominator. The
machine-readable decisions and exact scored atom identities are frozen in
[`benchmarks/atomic_reference_ledger.json`](../benchmarks/atomic_reference_ledger.json).
This audit does not replace the earlier site-only ledger or change any released
coordinates.

These are operational development-reference decisions. Passing means that the
specified modeled atoms satisfy the declared experimental-fit, correspondence,
chemistry and context rules. It does not establish a unique conformation,
independent atomic truth, biological specificity, complete biological assembly,
confirmatory performance or an unseen sequence/fold split. A partial-core score
cannot be reported as full-ligand atomic recovery.

## Timing and source separation

The auditor read structural provenance and source archives, and did not open
predicted poses, docking scores, pocket results or recovery metrics. The same
seven already released cases and the same ligand copies were retained; no crystal
copy was substituted for a better fit. The quantitative policy was written before
this auditor inspected any new density map and before any new pose metrics were
reviewed. The coordinating agent accepted the stated operational gates. Manual
source adjudications and exact scored atom subsets were then frozen before new
pose evaluation.

The policy SHA256 is
`93887a57a341a8b40eb1f755a4617570ee09fe703f6c4213ca891e5d5ef8f733`.
The final ledger SHA256 is
`8ed3887e2575f8ed92d73562f9b79c088036f326ed177af2c71b55522eedf383`.
The separately pinned manual adjudication SHA256 is
`8d6c8dfe5daf6718038c11944bcadd77db19d66449885cab9894de578d471fd1`.
A later source-driven amendment must be explicit and precede metric review; it
must not quietly overwrite these decisions.

## Evidence and operational gates

For all seven cases, the audit downloaded current deposited coordinates,
structure factors, wwPDB validation XML, and both weighted validation map
coefficient files. It reconstructed full-cell 2mFo−DFc and mFo−DFc Fourier maps
with Gemmi 0.7.5 at sampling rate 3, then subtracted the full-cell mean and divided
by its standard deviation. The [RCSB map documentation](https://www.rcsb.org/docs/general-help/electron-density-maps-and-coefficient-files)
describes these deposited validation coefficients and their use. Raw structure
factors are archived, but this audit did not independently re-refine the models
or produce ligand-omit maps. These are model-phased validation reconstructions;
model bias remains possible.

Each scored whole sugar residue must have:

- RSCC ≥ 0.80, RSR ≤ 0.25 and average occupancy ≥ 0.90 in the selected-copy
  validation record.
- At least 90% of its modeled heavy atoms with 2mFo−DFc density ≥ 1σ somewhere
  within a 0.50 Å sphere of the atom. The search uses the 33 points on a Cartesian
  0.25 Å lattice inside that sphere, including its center.
- Atom-center mFo−DFc density strictly greater than −3σ for every modeled atom.
- Exact released/deposited atom-name and coordinate correspondence, a consistent
  CCD/link graph and stereochemistry, and no unresolved chirality flag or
  absolute reported Mogul bond/angle Z-score ≥ 5.

These are declared conservative selection rules for this development audit,
not universal wwPDB quality standards, empirically calibrated uncertainty bounds
or a guarantee that every retained atom is individually resolved. The 90% rule
explicitly permits some weak atoms; the 0.5 Å search permits local displacement.
Reported geometry/torsion/clash flags below the exclusion gates remain in the
ledger. Rare torsions alone do not establish an invalid bound conformation.

The stereo check compares the signed volumes formed by consistently named CCD
heavy-atom neighbors against CCD ideal coordinates. At glycosidic centers, an
absent leaving oxygen is replaced by the deposited external acceptor oxygen in
that comparison. It detects handedness correspondence, rather than performing
a new refinement or automatically validating all carbohydrate conformations.
The [RCSB ligand-quality documentation](https://www.rcsb.org/docs/general-help/ligand-structure-quality-in-pdb-structures)
explains why experimental fit and chemical geometry require separate review.

A full-ligand endpoint requires all intended carbohydrate residues and all
released modeled heavy atoms to pass. A core is the largest covalently connected
component of passing **whole residues**, with at least one complete sugar ring.
The deterministic tie-break is heavy-atom count, then the ascending prepared
residue-number tuple. Individual weak atoms are never deleted to obtain a pass.
Manual unresolved link-geometry adjudications are separately identified.

The full parent ligand graph is preserved in subsequent preparation. Core atoms
are only a partial scoring mask: glycosidic leaving oxygens are not restored,
excluded residues are not removed, and chemical defects in the parent are not
silently repaired. In particular, 1QFO's accepted coordinate core does not make
its wrong-handed parent GLC chemically valid.

## Frozen decisions

Deposited author identifiers are shown below. Prepared ligand identifiers and
every exact heavy-atom name are in each case's `scored_atom_map`; each row is tied
to the coordinate archive SHA and released ligand compressed/uncompressed SHAs.

| Case | Qualified scope | Exact scored deposited residues | Atoms | Source evidence and retained limits |
|---|---|---|---:|---|
| 154L | Rejected | None | 0 / 43 | All NAG B1–3 fail RSCC (0.602/0.794/0.668) and 90% density coverage (12/15, 12/14, 11/14). B1 and B2 also have reported angle Z −10.53 and −5.32. EDS R 0.3099 differs substantially from depositor R 0.159. The published cleft remains a site-only reference. |
| 1GNY | Full modeled xylopentaose | XYP B1–5 | 46 / 46 | RSCC 0.808–0.971; occupancy 0.95–1. All whole residues pass. B1 is exactly at 9/10 supported atoms, with terminal O1 below the local 1σ rule. B5 has the lowest RSCC. Geometry flags up to absolute bond Z 4.68 remain; passing does not remove terminal uncertainty. |
| 1KJL | Modeled GAL core | GAL B2 | 11 / 26 | GAL has RSCC 0.902 and 11/11 density-supported atoms. NAG B1 has RSCC 0.751, only 13/15 supported atoms and O7 difference density −4.627σ. EDS R 0.3104 versus depositor R 0.186 remains an entry-wide warning even for the accepted local core. |
| 1QFO | Modeled sialyl-GAL core | GAL E2 + SIA E3 | 31 / 43 | Both pass local gates (RSCC 0.932/0.934; all atoms supported). GLC E1 has favorable density but wwPDB reports C1 WRONG HAND, independently reproduced by the CCD signed-volume check. Full ligand is rejected. Prepared receptor A is deposited protein B, and prepared ligand B is deposited glycan E. |
| 2J1T | Conditional modeled GAL–FUC core | GAL C2 + FUC C3 | 21 / 46 | Numerical fit, density and chirality gates pass all four residues. Manual unresolved link-geometry review excludes NAG C1 and FUC C4; the accepted core has RSCC 0.956/0.973 and all 21 atoms supported. No new native repair or invented bond Z-score is used. |
| 1PWB | Full modeled maltose | GLC D1–2 | 23 / 23 | RSCC 0.849/0.961 and all 23 atoms supported. D1 O1 difference density −2.867σ is close to the exclusion boundary; nearby negative difference features remain in the contour QA. Two clash flags and an angle flag remain. The released receptor is one chain of the experimental trimer. |
| 2VUZ | Modeled branched core | BMA B3 + MAN B4 + NAG B5 + MAN B6 | 47 / 90 | Residues B1/B2 have only 11/15 and 11/14 density-supported atoms; B1 RSCC is 0.700. B7 has angle Z −5.49. The largest remaining connected component passes, with B5 exactly 13/14 supported atoms. The administered nonasaccharide-Asn has only seven modeled sugars; two Gal and Asn are absent. |

For 2J1T, deposited `struct_conn` records `covale1` NAG C1 O4–GAL C2 C1
at 1.413 Å, `covale2` NAG C1 O3–FUC C4 C1 at **1.613 Å**, and
`covale3` GAL C2 O2–FUC C3 C1 at 1.418 Å. The unusually long second
C–O bond leaves native link geometry unresolved despite positive local density
and the paper's model assignment. Both whole residues touching that bond are
excluded. The within-residue validation flags for those residues are retained;
none supplies a measured interresidue ideal/sigma. The audit therefore assigns
no synthetic Z-score to that bond and does not pretend this manual adjudication
follows automatically from the numerical density gates.

## Primary publications and access scope

[Weaver et al. (154L)](https://doi.org/10.1016/S0022-2836(95)80038-7),
[Szabo et al. (1GNY)](https://doi.org/10.1074/JBC.M109558200) and
[Sörme et al. (1KJL)](https://doi.org/10.1021/ja043475p) were reviewed at
primary-abstract scope; a paper density figure was not independently viewed for
these cases. [Shrive et al. (1PWB)](https://doi.org/10.1016/S0022-2836(03)00761-7)
was accessible as publisher abstract/section excerpts. The new local-density
evidence for these cases comes from actual downloaded validation coefficients,
not an assumption that the abstract certifies exact atoms.

[May et al. (1QFO)](https://www.sciencedirect.com/science/article/pii/S1097276500800714)
reports clear complete sialyllactose density for molecules A/B; Fig. 1A illustrates
the already selected B-copy ligand at 1σ. Its C copy has only sialic acid modeled
and was not substituted. [Boraston et al. (2J1T)](https://www.jbc.org/article/S0021925819350045/pdf)
reports well-defined Lewis Y density and Fig. 5D at 1σ. These indexed publisher
full-text passages/captions were accessible, but direct raw full-text/PDF downloads
failed. [Gourdine et al. (2VUZ)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2662075/)
was reviewed in downloaded full text, including Methods, Results and Fig. 4;
it explicitly describes the unmodeled capping galactoses and Asn. The source
manifest distinguishes all access scopes and failures.

## Receptor and assembly context

Missing protein/cofactor atoms were checked against the actual selected ligand
and accepted subset. The 1PWB trimer's other deposited chains lie more than 22 Å
from selected maltose D, so their absence does not remove a direct atomic contact
at this site. Released chain A retains recognition calcium but selects one of
the two 0.5-occupancy Arg343 conformers; the omitted conformer also approaches
the ligand. This is a single-chain holo receptor experiment, not a complete SP-D
assembly claim.

2VUZ's publication describes a disulfide dimer while current assembly 1 is
annotated monomeric. Its Cys44-linked partner under operator `4_555`
(`-x,y,-z+1/2`) lies at least 15.696 Å from the retained core. The released
receptor omits HYP91; its nearest atom is 4.052 Å from the core but its carbonyl
oxygen lies 3.873 Å from excluded NAG B7. Ser61 alternate conformers also differ
near the core. These are retained receptor limitations, not an assertion that
a 4 Å cutoff removes all physical effects. Other constructs/domains and omitted
physiological context remain as recorded in the earlier structural provenance.

## Reproduction and archives

[`benchmarks/sources/atomic_reference_sources.json`](../benchmarks/sources/atomic_reference_sources.json)
records URLs, retrieval status, checksums, licenses and individual paper scopes.
All raw experimental sources, diagnostics, policy, manual decisions and seven
visually reviewed ring-plane contour figures are under
[`data/archive/atomic_reference/`](../data/archive/atomic_reference/).
Each contour is a slice through the sugar-ring mean plane with projected atom
labels; off-plane substituents must be assessed from the three-dimensional
per-atom values, not inferred from that slice alone.

The original source-only computations are retained as `inspect_sources.py`,
`inspect_receptor_context.py`, `render_density_slices.py` and `freeze_ledger.py`
in that archive. The reusable read-only replay is
[`scripts/verify_atomic_references.py`](../scripts/verify_atomic_references.py):

```sh
.venv/bin/python scripts/verify_atomic_references.py
.venv/bin/python -m pytest tests/test_atomic_reference.py
```

The actual read-only replay passed all seven cases, checking 86 archive hashes,
317 released source atoms, 26 sugar residues and 179 scored atoms; the maximum
normalized density difference was 0.0σ. Eight tests passed, including atom-loss
and moved-ligand controls whose source checksum pins were deliberately updated,
so those controls exercise inventory/coordinate replay beyond a checksum check.
The archived [`replay_verification.json`](../data/archive/atomic_reference/replay_verification.json)
is a derived output and is excluded from its own replay input inventory.

Raw-source checksums, atom inventories, identities, coordinates and cutoff
decisions are exact. Only recomputed density floats admit the explicitly stated
1e−5σ absolute tolerance, to accommodate Fourier/interpolation implementation
roundoff; this tolerance does not relax a qualification threshold. Separately
hashed manual judgments remain review decisions, not automated experimental
proof. Tampered raw maps, a lost ligand atom and moved ligand coordinates must
fail verification.

The [wwPDB usage policy](https://www.wwpdb.org/about/usage-policies) places PDB
archive data under CC0 1.0 and encourages citation of original structure authors.
Published articles retain their separately recorded publisher licenses. The
already archived upstream protein-glycans license is Apache 2.0; source metadata
and underlying structural provenance are retained.
