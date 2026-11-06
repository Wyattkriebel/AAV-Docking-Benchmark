# Source qualification of the twelve fixed candidate references

The source-only audit of 3 November 2026 retains all twelve candidates and
qualifies **six full modeled ligands, five declared modeled cores, and one atomic
rejection**. Eight candidates also meet the declared source-conditioned
preparation requirements; preparation remains blocked for 2G7C, 2GUP, 2WAB and
3A2S for the distinct reasons below. These structural decisions do not establish
family, fold or learned-model training independence. Those are separate audits.
The filename “heldout” identifies the candidate workflow, not an earned
independence claim.

The frozen ledger is
[`benchmarks/heldout_reference_ledger.json`](../benchmarks/heldout_reference_ledger.json),
SHA256 `5e1142c49191ee9881efc1e918630abc0f57cffb5f2fb0041350cd6df0ac7c4e`.
It contains exact standard CCD atom identities, formal charges, aromatic/stereo
metadata, the complete deposited parent bond graph, scored masks, biological
assembly transforms and every qualified copy. No new docking results, predicted
poses or recovery metrics were consulted. No old atomic-reference files or
frozen candidate sources were altered.

## Fixed policy and timing

The twelve-case denominator is the unchanged
[`heldout_candidate_inventory.json`](../benchmarks/heldout_candidate_inventory.json).
No case was replaced because of a method's outcome. The audit reuses the unchanged
operational policy at
`data/archive/atomic_reference/qualification_policy.json`, SHA256
`93887a57a341a8b40eb1f755a4617570ee09fe703f6c4213ca891e5d5ef8f733`.
All scored masks, original-reagent judgments, assembly/cofactor choices and
same-chemical-copy decisions were frozen before new runs.

Each retained whole sugar residue must satisfy RSCC ≥ 0.80, RSR ≤ 0.25, the
existing implementation's wwPDB residue-average occupancy (`avgoccu`) ≥ 0.90,
≥90% modeled heavy-atom density coverage, atom-center difference density >−3σ,
CCD/coordinate handedness correspondence, and no unresolved chirality flag or
reported absolute Mogul bond/angle Z ≥5. The prior JSON expresses an occupancy
minimum under residue-level requirements; its wording must not be overstated as
an explicit per-atom minimum. Individual occupancies and alternate rows remain
visible. Missing required diagnostics are not passes.

Coverage is measured from the actual weighted validation 2mFo−DFc map: at least
1σ at one of the 33 Cartesian 0.25 Å lattice points inside a 0.50 Å sphere about
each atom. Full-cell maps are reconstructed with Gemmi 0.7.5, Fourier sampling
rate 3, mean-subtracted and divided by the full-cell standard deviation. The
mFo−DFc atom-center gate is evaluated independently. These are operational
selection rules, not universal wwPDB standards or calibrated coordinate errors.

All twelve coordinate/validation sources, structure factors and both coefficient
files were available. Actual map sampling and visually reviewed ring-plane
contours supply local evidence; downloaded structure factors were not used to
independently re-refine models. The maps are model-phased and can retain model
bias. Published ligand-omit figures corroborate selected cases but are not newly
reconstructed omit maps. Passing does not prove every atom individually resolved,
a unique conformation, biological specificity or confirmatory performance.

Full qualification requires all deposited parent residues, including a covalent
aglycone, to pass and the original carbohydrate reagent to be complete. Otherwise
the core is the largest connected component of passing **whole sugar residues**,
with a complete sugar ring. Ties use heavy-atom count then lexical residue IDs.
An aglycone alone cannot form a glycan core. Individual bad atoms are never
removed to improve a pass, and alternative ligand conformers are not arbitrarily
selected. Core masks do not authorize truncating or repairing the preparation
parent.

## Twelve-case decisions

Residues below use current deposited label chains and branch residue numbers.
Original depositor author numbering and all standard atom names are separately
preserved in the ledger.

| Case | Structural scope and scored atoms | Retained source decision | Source-conditioned preparation |
|---|---|---|---|
| 1I3H | Full MAN B1–2, **23** | Both pass; retain the complete deposited tetramer, CA and MN. Four symmetry instances are one system. The known development fold overlap remains a separate generalization issue. | Eligible |
| 2G7C | Modeled NAG/GAL/GLA D1–3 core, **37** | All modeled sugars pass, but the primary abstract explicitly identifies an additional reducing-end O(CH₂)₈CO₂CH₃ spacer absent from the deposited graph. | Blocked: complete original parent unresolved |
| 2GUP | GLC B1 core, **11** | FRU B2 has difference density −3.488645σ. Source sucrose coordinates are not independent evidence that this structural-genomics site has the intended functional recognition interpretation. | Blocked: functional-site interpretation unresolved |
| 2RDK | Full MAN C1–2, **23** | Both pass. Assembly 1 is protein A/glycan C; source D belongs to distinct assembly 2. Retain the deposited construct state. | Eligible |
| 2WAB | BGC B1–3 core, **34** | B4 angle Z −5.30 fails; B5 has two 0.5-occupancy conformers. Only five of the original six sugars are modeled. | Blocked: missing original sugar and ambiguous parent alternatives |
| 2WNB | A2G B1 + GAL B2 core, **25** | Full parent is **35 atoms**, including covalently linked OPO author A402/label D. OPO has angle Z 6.12 and bond Z 8.32/−5.83, so its atoms are excluded from scoring. Preserve OPO in the chemical input and retain separate CMP (C5P) in the primary receptor state. | Eligible for the declared core |
| 3A2S | Rejected, **0** | Required RSCC/RSR/avgoccu are absent. Only 5/11 GLC and 7/12 FRU atoms meet density coverage; GLC also has −3.133165σ atom-center difference density. Experimental files exist; this is not an unavailable-map rejection. | Blocked: required atomic evidence fails |
| 3ZWE | Full BGC/GAL/FUC/GLA F1–4, **44** | F passes. Same-chemical E fails GLA4 RSCC and cannot enter the full 44-atom minimum. Complete protein A/B/C trimer retained; partial sites remain known-site annotations. Exact external training overlap is a separate restriction. | Eligible structurally |
| 4R9F | BMA B2 core, **11** | B1 O1 has −4.278822σ difference density despite good residue fit and published ligand-omit support. The fixed gate excludes the entire B1 residue. Whole two-domain monomer and explicit MSE identity retained. | Eligible for the declared core |
| 5HZA | Full BGC/FUC/GAL C1–3, **33** | All three 3-fucosyllactose residues pass; retain both P-domain chains A/B. This does not establish intact-capsid accessibility. | Eligible |
| 5V6F | Full BMA/MAN/MAN B1–3, **34** | All pass. Published mannotriose is predominantly one orientation; mixed glycan/glycerol occupancy discussed for another deposition is not assigned to this case. Sulfate is an inorganic additive. | Eligible |
| 8I5U | Full BGC B1–2, **23** | Both laminaribiose residues pass. This is the deposited mutant's noncovalent enzyme–substrate complex, not a protein–glycan covalent intermediate. Whole source monomer and CA retained. | Eligible |

The 4R9F inventory's expansion of “Cp” to *Clostridium perfringens* is incorrect;
the primary paper and deposited source identify *Caldanaerobius polysaccharolyticus*.
That correction is recorded in the new audit without modifying the fixed
inventory. The 2WAB source CIF explicitly identifies a mutant construct; a
wild-type density figure in its paper is not relabeled as the selected mutant's
experimental evidence.

## All-copy truth and receptor frame

The primary selected-copy mask is fixed first. Every deposited copy of the same
complete chemical graph is then independently checked, including CCD residue/link
identity and covalent aglycones. A copy enters atomic scoring only if it belongs
to the selected biological assembly and passes **the same canonical atom mask**.
There is no minimum over a full ligand in one copy and a smaller convenient core
in another.

The audit replays 19 same-chemical asymmetric-unit copies and freezes 16 qualified
assembly instances:

| Case | Qualified instances in selected assembly | Other source sites retained |
|---|---:|---|
| 1I3H | B under four deposited tetramer operators | Copies are symmetry-related, not independent observations. |
| 2G7C | C and D, two modeled-core instances | E/F are in assembly 2, outside the selected receptor. |
| 2GUP | B under two dimer operators | C/D fail the fixed GLC mask; D's passing FRU cannot substitute for it. |
| 2RDK | C, one | D is quality-passing but belongs to assembly 2 and is not unioned into assembly 1. |
| 2WAB | B, one core instance | Unscored B5 alternatives remain separate; neither is selected for parent preparation. |
| 2WNB | B plus attached D/OPO parent, one sugar-core instance | CMP remains a separate co-ligand. |
| 3A2S | None | Published sucrose site remains known, without atomic qualification. |
| 3ZWE | F, one | E fails the fixed full mask; D is a disaccharide and G/H/I are isolated fucose sites. |
| 4R9F, 5HZA, 5V6F, 8I5U | One each | Exact source ligand and assembly retained. |

`qualified_copies` contains literal deposited operator matrices/vectors and
transformed scored coordinates in the same world frame as the complete receptor.
Rigidly generated copies inherit the source copy's experimental quality; they are
not sampled against an untransformed crystal map. Atomic evaluation may minimize
over these exact copies after inputs and ranks are frozen. Binding at another
known but partial/quality-failing site is an **unscored known site**; it cannot be
called biological nonbinding merely because it lacks qualified atomic truth.
Neither selected sites nor these annotations may define prediction crops, boxes,
pockets or rankings.

The primary receptor policy keeps every protein subchain/operator in the selected
deposited construct assembly, native covalent modifications, all deposited CA/MN,
and 2WNB CMP. It removes waters, other glycans and declared crystallization
additives (GOL, EDO, sulfate, iodide, TRS). Protein alternatives follow a generic
highest-occupancy/lexical tie rule with exact source mapping. MSE is preserved as
MSE; no silent MSE→MET substitution is authorized. Source CCDs for C5P, MSE, CA
and MN are archived. Method feasibility and any proposed alternative preparation
must be explicit rather than redefining structural eligibility.

## Primary sources and access limits

The new [primary-paper evidence record](../data/archive/heldout_reference/papers/primary_paper_evidence.json)
provides individual URLs, reviewed scope, ligand identity and license context.
Full primary text was archived for the six priority cases:
[2RDK](https://pmc.ncbi.nlm.nih.gov/articles/PMC2327273/),
[3ZWE](https://pmc.ncbi.nlm.nih.gov/articles/PMC3281725/),
[4R9F](https://pmc.ncbi.nlm.nih.gov/articles/PMC4263893/),
[5HZA](https://pmc.ncbi.nlm.nih.gov/articles/PMC4836343/),
[5V6F](https://pmc.ncbi.nlm.nih.gov/articles/PMC5825169/), and
[8I5U](https://pmc.ncbi.nlm.nih.gov/articles/PMC10167657/).
Original density figures for 3ZWE, 4R9F, 5HZA and 5V6F were also archived and
inspected. These are source-only corroboration, not outcome-based selection.

2WAB and 3A2S primary full texts were accessible. 1I3H and 2G7C were reviewed at
primary-abstract scope; the 2G7C abstract itself states the missing spacer formula.
2WNB's official supplementary PDF identifies the nitrophenyl acceptor and separate
CMP; its main article was not obtained. 2GUP has deposition metadata without a
primary journal assignment. Failed PDF requests and HTML responses are honestly
recorded; an HTML response is not called a downloaded supplement.

Published omit maps do not waive the frozen current-source gates. For example,
the 3A2S paper supports a sucrose transport site and the 4R9F paper illustrates
ligand-omit density. The new decisions concern the exact current selected atom
sets under the declared policy, not whether those authors observed binding.

## Reproduction

The source manifest is
[`data/archive/heldout_reference/source_manifest.json`](../data/archive/heldout_reference/source_manifest.json).
It preserves the 92 original candidate pins, all new experimental retrievals,
raw/source checksums and separate article license contexts. Derived density PNGs
and replay output are recorded separately to avoid a report hashing itself.
The [wwPDB policy](https://www.wwpdb.org/about/usage-policies) places PDB archive
data under CC0; articles keep their individual licenses.

The public package retains **163 byte-pinned records**. Eleven publication or
provider-page snapshots remain metadata-only records with their original hashes,
byte counts and retrieval locations; permission to redistribute their bytes was
not established. Local research copies remain unchanged. Every coordinate,
validation, coefficient-map and chemical-definition input used numerically is
retained. The separate [public distribution manifest](../benchmarks/source_public_distribution.json)
records the exact omissions and attribution for the retained licensed articles
and figures. See the [independent packaging review](SOURCE_COHORT_PR_REVIEW.md).

```sh
.venv/bin/python scripts/verify_public_source_cohort.py
.venv/bin/python -m pytest tests/test_heldout_reference.py -q
```

The original local replay verified 174 source/archive hashes. Public replay
verifies all 163 retained source records, acknowledges the eleven noncomputational
publication records separately, and replays all twelve cases, **388
native heavy-atom rows representing 377 unique parent atoms**, 298 primary scored
atoms, 19 same-chemical source copies and 16 qualified assembly instances. The
11 extra rows are 2WAB's second terminal alternate; its 67 rows must not be
mistaken for 67 unique atoms or a fully resolved hexasaccharide. The largest
normalized density replay discrepancy is 0.0σ.

Fourteen tests passed, including native atom loss/movement after raw checksum
pins were deliberately updated, corrupted map bytes, lost scored-mask atoms,
misplaced transformed copies and substitution of a good copy from the wrong
assembly. Only normalized map floats admit the stated absolute 1e−5σ tolerance;
identities, inventories, coordinates, graphs, cutoff classifications, canonical
masks and assembly transforms are exact. Manual paper/reagent/context judgments
remain separately hashed decisions rather than claims proved by a numerical
replay.
