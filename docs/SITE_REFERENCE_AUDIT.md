# Published-site reference audit

Audit date: 30 October 2026. The eight existing development systems qualify for
**site-only localization of the published modeled glycan footprint**. None is
promoted to an independently certified atomic-pose reference, specificity panel,
confirmatory test case, or independent fold split. Decisions and case cautions
were communicated before the auditor inspected any new pocket outcomes; no new
pocket outcomes were consulted. The cohort is AAV9 plus every case already frozen
in `benchmarks/external_cohort.json`, without selecting cases for a method's result.
Machine-readable decisions and evaluation-only identifiers are in
`benchmarks/site_reference_ledger.json`; primary citations and review scope are in
`benchmarks/sources/site_reference_sources.json`.

## Criteria fixed before outcome review

Site-only acceptance requires all of the following:

1. A primary structural publication assigns a noncovalent glycan-recognition
   surface or anchor, rather than merely reporting a carbohydrate somewhere in a
   coordinate file.
2. The modeled deposited ligand corresponds to the pinned reference, with the
   selected chain and residue identities explicit.
3. Deposited validation is reviewed for that ligand copy, retaining local-fit,
   occupancy, geometry and entry-wide warnings. Missing evidence stays missing.
4. The receptor representation and relevant omitted construct/subunit context are
   disclosed. Site prediction uses the whole pinned released bound receptor, or
   the whole deposited AAV capsid, without ligand-defined crops or restrictions.
5. The endpoint is approximate site localization. Acceptance does not certify
   unseen local density, a chemically valid exact atomic pose, biological
   specificity, full biological assembly, or generalization to unseen folds.

There is no universal resolution, RSCC or occupancy cutoff. Decisions use each
paper's site assignment, coordinate correspondence, and the selected deposition's
uncertainties. Favorable RSCC does not override chirality, geometry or assembly
problems; poor exact-atom support need not erase a well-supported recognition
surface. No local electron-density or cryo-EM maps were independently inspected.
The original atomic-reference gates in `PROTOCOL.md` and the earlier reference
ledger remain in force.

The existing protocol's center-to-nearest-native-heavy-atom Top-1/Top-5 endpoint
uses 4 Å with 3/5 Å sensitivity, together with modeled-interface residue coverage.
Those distances define the operational localization comparison; they are not
confidence bounds on the experimental atom positions. For incomplete glycans,
this evaluates only the modeled footprint. A site-only result cannot be
reinterpreted as accurate full-ligand docking.

## Case decisions and retained uncertainty

The numerical diagnostics below come from the already archived
`benchmarks/sources/reference_validation.json`. They describe the **selected
reference ligand**, not whichever symmetry or crystal copy has the best metric.
The linked sources are the original structural publications. Where only their
abstract was accessible/reviewed, the source record says so.

| Case | Why the published site is usable | Why atomic/context certification remains unresolved |
|---|---|---|
| AAV9 / 7MTZ | [Penzes et al. 2021, Results/Fig. 4](https://pmc.ncbi.nlm.nih.gov/articles/PMC8428384/) locate ligand-dependent density in the recognition pocket near a threefold spike. All 60 modeled GAL sites belong to one equivalent site class. | The paper calls glycan density poorly defined. Methods dose GalNAc while the deposited component is GAL. The article associates 7MTZ with 2.7 Å; current deposition says 2.43 Å. Archived validation provides no local ligand-fit metric and flags one torsion per GAL copy. No exact atomic pose is certified. |
| 154L | [Weaver et al. 1995](https://doi.org/10.1016/S0022-2836(95)80038-7) place noncovalent chitotriose in the B/C/D catalytic subsites. All 43 reference ligand atoms match deposited coordinates exactly. | RSCC 0.602–0.794, occupancy 1, 12 angle, four torsion and one bond flag; EDS R 0.3099 versus depositor R 0.159. This is a lower-confidence modeled atomic footprint within a published catalytic cleft. |
| 1GNY | [Szabo et al. 2001](https://doi.org/10.1074/jbc.M109558200) identify a 20–25 Å concave groove binding xylopentaose. All 46 modeled atoms map exactly. | RSCC 0.808–0.971 and occupancy 0.95–1; eight angle and five bond flags. Terminal support is weaker. The isolated CBM15 is not the complete Xyn10C enzyme. |
| 1KJL | [Sorme et al. 2005](https://doi.org/10.1021/ja043475p) solve the CRD–LacNAc complex. The GAL anchor (RSCC 0.902) localizes its recognition surface; all 26 ligand atoms map exactly. | NAG RSCC 0.751, two clash and one bond flag; EDS R 0.3104 versus depositor R 0.186. One prepared terminal oxygen differs from deposition, outside the site. This CRD omits the full-length protein's N-terminal region. |
| 1QFO | [May et al. 1998](https://doi.org/10.1016/S1097-2765(00)80071-4) assign the functional sialic-acid-binding template. Exact mapping selects protein author B and ligand author/label E, whose three residues have RSCC 0.928–0.934 and occupancy 1. | Selected GLC has a chirality flag and selected SIA has two angle flags. Weaker isolated SIA C201 is a different reference. This is only the N-terminal receptor domain. Site use does not repair the full ligand's stereochemistry. |
| 2J1T | [Boraston et al. 2006](https://doi.org/10.1074/jbc.M607620200) describe the CBM47 Lewis Y recognition surface. All 46 atoms map; selected ligand C1–4 has RSCC 0.899–0.973, occupancy 1. | Nine angle and two torsion flags remain in the four selected residues. The exact flags/local density have not been independently certified. One module does not represent the complete multivalent parent protein. |
| 1PWB | [Shrive et al. 2003](https://doi.org/10.1016/S0022-2836(03)00761-7) establish maltose recognition by recombinant SP-D. The selected D1–2 ligand maps exactly; RSCC 0.849–0.961, occupancy 1. Recognition calcium remains in the released source. | Two clash and one angle flag in selected maltose. The paper and deposited assembly 1 describe a trimer; the released receptor has one chain. It cannot support a complete-assembly claim. |
| 2VUZ | [Gourdine et al. 2008](https://pmc.ncbi.nlm.nih.gov/articles/PMC2662075/) identify calcium-mediated glycan recognition and explicitly report unmodeled capping galactosides. All 90 atoms of the seven modeled sugars map exactly. | The title's nonasaccharide has only seven modeled residues. RSCC 0.700–0.939, occupancy 1, 18 angle, five torsion and two bond flags. Released protein omits HYP91. Published dimer context conflicts with current monomer assembly annotation. |

The 2J1T angle total is **nine** when summing the four selected validation records
in the archived JSON; the earlier `NEXT_EVIDENCE.md` prose says eleven. The audit
uses the explicit selected records, without silently carrying over that prose
count. Likewise, the 1QFO selected-copy assessment does not pool the isolated
sialic acid or the alternative crystal copy. Whole-entry warnings are retained
separately and are never removed by copy selection.

## Receptor and assembly scope

“Full released bound receptor” means every atom in the pinned upstream receptor
file, without a native-site crop. It does not mean full-length physiological
protein or independently complete biological assembly. For 154L, 1GNY, 1KJL,
1QFO and 2J1T, the released protein chain corresponds to the selected deposited
monomeric construct (1QFO assembly 2; otherwise assembly 1). Alternate conformers,
waters and some source preparation details still differ. A holo structure can
carry binding-site conformational information even after its ligand is removed.

1PWB is explicitly a single-chain representation of a trimeric experimental
construct. 2VUZ's publication describes a disulfide-linked dimer; its current
mmCIF assembly 1 is author/PQS annotated monomeric, while `struct_conn` records a
Cys44–Cys44 disulfide to crystal operator `4_555`. That discrepancy is recorded,
not resolved by assuming either annotation establishes the entire biological
context. Its released receptor also drops HYP91 and renumbers later residues.

For AAV9, the deposited model already contains 60 protein chains and 60 GAL
molecules; assembly 1 uses the identity operator over that expanded model. The
prediction receptor must contain the full deposited capsid exactly once. This
is the complete **deposited ordered capsid**, with disordered VP regions and the
full host glycoconjugate absent. Sixty equivalent sites are one observation, not
60 independent successes.

## Evaluation identities and prevention of native-site leakage

The ledger's `evaluation_only` fields preserve exact author/label ligand and
contact identifiers. They are for scoring only after ligand-free prediction
inputs and candidate ranks have been frozen. Native coordinates, residue lists,
mutagenesis sites, papers and validation information must not select a prediction
crop, sector, pocket, box, surface patch or random-baseline distribution.

All seven external ligand files exactly match the deposited atom names and
coordinates at 0.001 Å precision. This verifies correspondence, not chemistry or
density. Protein contact identities are mapped using exact atom coordinates
rather than a guessed uniform numbering offset. In particular, released 1QFO A
is deposited protein B and released ligand B is deposited glycan E. 2VUZ needs a
nonuniform map because HYP91 was omitted. The source's renamed calcium atoms are
reported separately from protein-residue contacts.

Contacts use the protocol's geometric rule: any protein heavy atom strictly
within 4 Å of a modeled native ligand heavy atom. This differs from a paper's
biochemical description of a direct interaction. For example, GAL author A801
(label IB) in 7MTZ geometrically contacts Trp503 on author A and Asn470/Val473 on
neighboring author I, while Penzes et al. describe Trp503 as the direct glycan
contact. Both statements are retained with their definitions. Monomer evaluation
would omit the neighboring-chain contacts.

## Source families and allowed conclusions

The source labels describe goose-type lysozyme, CBM15, galectin-3, a siglec V-set
domain, CBM47/fucolectin-related protein, SP-D, codakine and AAV capsid. These are
useful descriptions of this small development collection. They are not an audit
of independent sequences or folds: SP-D and codakine are both C-type lectins,
and several other entries have broad beta-sandwich architectures. Homolog and
training-set overlap remain unaudited. The existing external-cohort label calling
1GNY an endo-xylanase family must not obscure that the actual released receptor is
its noncatalytic CBM15 domain.

The earned result is a transparent, fixed eight-system **development comparison
of published-site localization**, with model-derived contact coverage and the
listed quality/context caveats. Atomic RMSD can remain an implementation or
source-reproduction diagnostic, but this audit does not upgrade it to validated
atomic accuracy. No new specificity, complete-assembly or held-out generalization
claim follows from the site-only decisions.
