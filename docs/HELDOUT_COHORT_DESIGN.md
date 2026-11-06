# Source-frozen held-out cohort candidate audit

The [fixed inventory](../benchmarks/heldout_candidate_inventory.json) records the
initial source-availability snapshot of **12 candidate constructs**. Its frozen
zero-acceptance fields describe that initial stage. The subsequent
[atomic qualification](HELDOUT_REFERENCE_QUALIFICATION.md) accepts six full
modeled ligands and five declared modeled cores, and rejects one atomic
reference. The completed [family/training audit](FAMILY_LEAKAGE_AUDIT.md)
certifies **zero independent novel or learned-method-held-out families**.
The inventory was selected from source
coordinates, chemistry, validation reports and released benchmark metadata on
3 November 2026 before new docking for these cases in this repository. It is an
availability audit that makes the next qualification work concrete; it does not
close the confirmatory study in [PROTOCOL.md](../PROTOCOL.md).

The fixed six-case starting set is **2RDK, 3ZWE, 4R9F, 5HZA, 5V6F and 8I5U**.
All six now have source-qualified atomic endpoints: five full modeled ligands
and the 11-atom BMA B2 core of 4R9F. The original inventory separately lists
seven priority map-review candidates, including 1I3H; its proposed six-case
starting set excludes that known development-superfamily overlap. Neither list
was replaced in response to predictions. Local maps, original reagents,
assemblies and same-chemical source copies have now been audited. Family and
training exposure remain explicit restrictions on what these cases can test.

## Fixed sources and selection boundary

The released [HADDOCK-ready benchmark](https://github.com/haddocking/glycan-docking-benchmark/tree/765ac7728869ef47611c6aa68f53f2adcfa39207/HADDOCK-ready)
provides bound/unbound inputs and original reference structures. The source tree
is pinned at `765ac7728869ef47611c6aa68f53f2adcfa39207`. This audit uses its
availability and ligand descriptions, not its predictions or precomputed RMSDs.
The original development inventory remains unchanged.

[ProtGlycanDock](https://arxiv.org/html/2610.26009v1) and its
[released dataset](https://huggingface.co/datasets/ProtGlycanDock/ProtGlycanDock/tree/fbbb866c7d0909391364d93faa85b87b7c5b6db7)
provide additional coordinate/checksum and split metadata. The authors' protein
holdout is not automatically a holdout from our development cases, all structural
folds, P2Rank, GNINA or other pretrained models. A page retrieval unexpectedly
included published docking outcome excerpts; the named demonstration case 2JEN
was excluded. Public benchmark reuse is not a claim of universal blinding. No
old repository docking outcome files were used to choose this inventory.

Every candidate remains fixed in PDB-ID order. Source-based exclusions and their
reasons must be locked before predictions. An unsuitable candidate is recorded
as such; it must not be silently replaced after outcomes are seen. Subsequent
new candidates require a separately dated source-only inventory.

## Initial availability snapshot and its subsequent qualification

The table below preserves the availability questions recorded before the later
map and family audits. It is not the current acceptance ledger. Current masks,
rejections, source-conditioned preparation eligibility and all qualified copies
are in [the qualification report](HELDOUT_REFERENCE_QUALIFICATION.md); completed
sequence/domain comparisons are in [the overlap audit](FAMILY_LEAKAGE_AUDIT.md).

The table uses current `_atom_site` author chain/residue IDs, which can differ
from the original depositor IDs in `_pdbx_branch_scheme`. Both maps, the selected
branch entity, every expected heavy-atom name, exact glycosidic endpoints,
validation values, source revisions and coordinate-identity hashes are retained
in the JSON. RSCC ranges summarize the selected sugars only. The later atomic
review applies its separately frozen operational RSCC and density gates; these
are not asserted to be universal structural-validation standards.

| PDB and proposed family | Selected deposited graph; heavy atoms | Resolution; sugar RSCC | Qualification decision still needed |
|---|---|---|---|
| [1I3H](https://www.rcsb.org/structure/1I3H), legume lectin/ConA | B1–2, α1→2 dimannose; 23 | 1.20 Å; 0.916–0.952 | Exact CATH `2.60.120.200` overlap with development galectin 1KJL; cannot count as a strict new fold. Retain Ca/Mn and review tetramer. |
| [2G7C](https://www.rcsb.org/structure/2G7C), TcdA CROPS | D1–3, Galα1→3Galβ1→4GlcNAc; 37 | 2.00 Å; 0.939–0.974 | [Primary source](https://pubmed.ncbi.nlm.nih.gov/16622409/) specifies `O(CH2)8CO2CH3` reducing-end spacer absent from deposited graph. Sugar completeness is not complete original ligand identity. |
| [2GUP](https://www.rcsb.org/structure/2GUP), ROK family | B1–2, sucrose; 23 | 2.01 Å; 0.873–0.892 | No journal citation in current entry; three sucrose copies, differing occupancy. Establish biological recognition versus additive/crystal contact. Explicit SeMet preparation map needed. |
| [2RDK](https://www.rcsb.org/structure/2RDK), cyanovirin-N | C1–2, α1→2 dimannose; 23 | 1.35 Å; 0.960–0.969 | Priority. Five-site mutant and tag; adjudicate selected monomer/ligand assembly and map. Inference is to that construct. |
| [2WAB](https://www.rcsb.org/structure/2WAB), carbohydrate esterase 2 | B1–5, β1→4 glucopentaose; 56 | 1.90 Å; 0.822–0.981 | Titled hexasaccharide is incomplete; terminal sugar has 0.5/0.5 A/B alternates. N-terminal CATH `2.60.120.260` overlaps development CBM15/CBM47-like domains. |
| [2WNB](https://www.rcsb.org/structure/2WNB), GT29 sialyltransferase | B1–2 plus A402 OPO aglycone; 35 | 1.55 Å; 0.950–0.965 | Full Galβ1→3GalNAc–O–ortho-nitrophenyl graph is deposited. A bare disaccharide would discard the covalent aglycone. Freeze CMP/cofactor and SeMet policies. |
| [3A2S](https://www.rcsb.org/structure/3A2S), PorB porin | A1–2, sucrose; 23 | 2.20 Å; absent | Validation XML has no local sugar RSCC/RSR values. Obtain usable maps and review membrane trimer, detergent and binding-site relevance. |
| [3ZWE](https://www.rcsb.org/structure/3ZWE), BambL/RSL propeller | F1–4, blood-group-B tetrasaccharide; 44 | 1.75 Å; 0.905–0.963 | Priority. Whole trimer supplies recognition context. Other branch E has terminal RSCC 0.720; branch D is a fragment. Exact PDB occurs in ProtGlycanDock training. |
| [4R9F](https://www.rcsb.org/structure/4R9F), CpMnBP1 SBP-II | B1–2, β1→4 mannobiose; 23 | 1.40 Å; 0.951–0.974 | Priority. Tryptic, SeMet-labelled construct; keep both binding domains and document atom changes. |
| [5HZA](https://www.rcsb.org/structure/5HZA), norovirus P domain | C1–3, 3-fucosyllactose; 33 | 1.35 Å; 0.925–0.988 | Priority. Keep complete P-domain dimer; inference to isolated P domain does not establish complete-capsid accessibility. |
| [5V6F](https://www.rcsb.org/structure/5V6F), RbmC β-prism | B1–3, branched α1→3/α1→6 mannotriose; 34 | 1.50 Å; 0.945–0.980 | Priority. Isolated second β-prism domain, not intact RbmC. Inorganic SO4 is an additive, not a sulfated glycan. |
| [8I5U](https://www.rcsb.org/structure/8I5U), GH116 | B1–2, β1→3 laminaribiose; 23 | 1.40 Å; 0.958–0.968 | Priority. D593N acid/base mutant; assign recognition domain and compare its other domains to development folds. Largest receptor, 799 canonical residues. |

All selected deposited graphs match CCD expected heavy atoms after explicit
leaving-atom removal, including the OPO linkage in 2WNB. This check establishes
atom-set consistency. It does not establish original-reagent completeness,
anomer/stereochemical correctness, acceptable ring/linkage geometry, absence of
model bias or interpretable density. In particular, it cannot rescue 2G7C or
2WAB as whole original-ligand references.

Primary articles are linked per case in the inventory. Additional useful source
descriptions include the [BambL assembled recognition sites](https://pmc.ncbi.nlm.nih.gov/articles/PMC3281725/),
[CpMnBP1 construct and ligand density](https://pmc.ncbi.nlm.nih.gov/articles/PMC4263893/),
[RbmC domain structure](https://doi.org/10.1371/journal.ppat.1006841) and
[norovirus milk-oligosaccharide complexes](https://doi.org/10.1128/JVI.03223-15).

## Original source-qualification plan, now executed

1. Review the selected ligand and neighboring recognition atoms against local
   experimental maps, with a second independent adjudicator. Fetch and archive
   the deposited structure factors, map coefficients/maps and their hashes;
   coordinate and validation-report availability alone does not certify a map.
   Check occupancy, alternate conformers, anomer, stereochemistry, ring/linkage
   outliers and unresolved tails. The inventory's unviewed-map fields describe
   the initial snapshot; the subsequent ledger and twelve reviewed contour
   figures record the completed source-only map audit. No docking pose may
   determine the scored core.
2. Reconcile the complete original ligand with the primary paper and deposited
   connectivity. Preserve exact reducing-end modifications, noncarbohydrate
   aglycones, metal/water requirements and formal charge. Do not substitute a
   chemically convenient ligand or reconstruct missing native tails.
3. Materialize the declared recognition assembly once. Match every selected
   native atom to its assembly copy and archive the transforms. Establish whether
   the isolated construct supports the intended estimand. Keep source-based
   receptor preparation separate from evaluation-only ligand/site coordinates.
4. Freeze eligible families and cases, atom maps, receptor inputs, chemical
   graphs, methods, configurations, failures, estimands and all information
   regimes in a new qualification/run manifest. Only then begin new predictions.
   Archive any rejected candidates and reasons without outcome-driven replacement.

The six-case starting set has passed source qualification, retaining full versus
core endpoints separately. Eight of the twelve candidates meet the declared
source-conditioned preparation requirements; this does not guarantee software
support or successful method preparation. The fixed six support a small
prospective descriptive workflow check with the exposure strata disclosed, not
the protocol's 12-family confirmatory cohort. Structural eligibility was not
relaxed to obtain that number.

## Recognition domains and training overlap

Exact canonical deposited sequences, lengths, UniProt mappings and source
metadata hashes are archived for all candidates. The inventory also records
source-only domain annotations for development 154L, 1GNY, 1KJL, 1QFO, 2J1T,
1PWB, 2VUZ and AAV9. A new PDB ID or a different protein name is insufficient.

The completed screen compares whole chains and source-defined CATH/SCOP domains
against development and the official default P2Rank CHEN11-fpocket protein list.
Its separate [policy](../data/archive/family_audit/qualification_policy.json),
frozen before similarity calculation, uses at least 30% identity among aligned
nongap residue pairs and 70% coverage of the shorter sequence. The inventory's
earlier 30%/80% suggestion remains an unexecuted proposal in that historical
file; it is not the threshold used in the published audit. All 56,366 alignments,
warnings and source-defined boundaries are retained. Exact CATH superfamily
evidence remains separate from numerical warnings, including matches below the
sequence threshold. Neither a screen-negative result nor a missing annotation
certifies independence.

CATH topology/superfamily, SCOP2 and primary structural evidence must be audited
alongside sequence. Source annotations already show exact fold overlap for ConA
versus development galectin, and the CE2 N-terminal domain versus development
CBM15/fucolectin-like domains. GH116 has a six-hairpin catalytic domain and a
second galactose-mutarotase-like domain, so the actual recognition domain must
be resolved before claiming fold independence. Missing or discordant annotations
require structural searches/alignments. A strict fold-disjoint analysis and a
new-sequence-family analysis with shared folds are different strata.

The completed audit finds two development CATH-superfamily overlaps and seven
default P2Rank training CATH-superfamily overlaps, despite zero exact PDB-ID
matches to that training list. The fixed six therefore are not six
P2Rank-training-held-out families. Gaps in classifications and warning
interpretation remain unknown rather than evidence of novelty.

The pinned ProtGlycanDock manifests mark 3ZWE as **train**, 1I3H as **valid**, and
2GUP/3A2S/2WNB/8I5U as **protein_holdout_test**. Its 3ZWE row uses a different
branch; the protein itself is still exposed in training. Absence of an exact
PDB from those manifests does not prove family-level or historical pretraining
absence. The P2Rank protein-list audit is now complete; its feature construction
and score calibration, GNINA and other learned methods still require separate
provenance audits. Unknown model provenance remains unknown. Classical methods can be tested on
new repository targets without claiming learned-model training independence.

## Estimand, uncertainty and a bounded stopping point

After qualification, a defensible primary comparison is the paired difference
in family-macro Top-1 heavy-atom RMSD ≤2 Å success between two prespecified,
equally informed methods. Average seeds within complexes, then complexes within
families. Keep eligible method failures in the denominator. Report per-family
and per-complex outcomes, continuous RMSD, Top-5, geometry and runtime/failures.
Site discovery, conditional predicted-site docking and oracle-local docking
remain separate estimands; an oracle result is not blind-site accuracy.

Use the protocol's 10,000 whole-family paired bootstrap resamples and 95%
intervals, with the small number of families stated prominently. Bootstrap
intervals can degenerate at all-success/all-failure and do not certify precision.
An exact binomial interval is applicable only to a genuinely independent,
one-case-per-family binary endpoint, not to a fractional family-macro mean.
Even 12/12 independent successes have a two-sided exact 95% lower bound of about
0.735; 6/6 have a lower bound of about 0.541. These are illustrations of sampling
uncertainty, not performance targets, and do not account for reference,
selection or training leakage uncertainty.

Finish the fixed eligible cohort and all frozen conditions regardless of the
observed ranking or interval width, then stop with descriptive estimates and
limitations. Do not invent a confirmatory precision floor. Wider claims need an
independently selected enlargement or a new untouched cohort after any parameter
change. Twelve candidate labels do not supply twelve independent families.

For a bounded initial CPU plan, one bound construct per family and three frozen
seeds can precede any large factorial expansion. A prospective cap of 20 elapsed
minutes and two CPUs per Vina invocation would bound two methods ×12 cases
×3 seeds at 48 CPU-hours for one information regime. This is a planning cap,
not an achieved runtime or equality of algorithmic search effort. Validate
method settings/caps on development first and retain test timeouts as failures.
Apo, water, assembly and five-pocket arms need their own fixed budget manifest;
seeds and sites do not create additional independent families.

## Portable source archive and verification

The later qualification adds 174 pinned source records. Public distribution
retains 163 and records eleven noncomputational publication/context snapshots
as metadata only; the original local manifest, ledger and source bytes remain
unchanged. The [distribution manifest](../benchmarks/source_public_distribution.json)
identifies every omission by original hash, size, URL and reason and supplies
attribution/license records for retained open-licensed papers and figures.
Free public article access is not treated as a redistribution license.

The public verifier replays all numerical maps, graphs, masks, assembly copies
and family alignments without downloading those publications. Paper-based
manual judgments remain provenance for independent source review, rather than
claims re-proven by numerical replay:

```sh
python scripts/verify_public_source_cohort.py
```

The [source manifest](../data/archive/heldout_candidates/source_manifest.json)
pins 92 source byte records: 12 coordinate mmCIFs, 12 compressed validation
reports, 12 entry metadata snapshots, candidate sequence/domain metadata, 11 CCD components, development
sequence/domain metadata, the released dataset file manifest and all four split
manifests. The archive is about 2.6 MB. Non-gzip sources use deterministic gzip
with mtime zero; validation reports retain the exact downloaded gzip bytes.
Coordinate sources are current revised PDB entries, and revision histories are
recorded. Future download changes must fail checksum comparison and trigger
review, not silently update an accepted reference.

From the repository root, verify every archived source without network access:

```sh
python - <<'PY'
import gzip, hashlib, json
from pathlib import Path
manifest = json.loads(Path('data/archive/heldout_candidates/source_manifest.json').read_text())
for item in manifest['sources']:
    archived = Path(item['archive_path']).read_bytes()
    assert hashlib.sha256(archived).hexdigest() == item['archive_sha256'], item['filename']
    source = archived if item['filename'].endswith('.gz') else gzip.decompress(archived)
    assert len(source) == item['bytes'], item['filename']
    assert hashlib.sha256(source).hexdigest() == item['sha256'], item['filename']
print(f"Verified {len(manifest['sources'])} pinned source byte records")
PY
```

Graph checks used Gemmi 0.7.5 to parse the deposited mmCIF/CCD, matching heavy
atom identities by label chain, current author residue number and CCD atom name,
then removing only declared heavy leaving atoms. The special 2WNB OPO–A2G C1
connection is retained explicitly. Full CCD source atom/bond/stereochemistry
records remain archived for independent re-adjudication. No new docking was run
as part of this availability audit.
