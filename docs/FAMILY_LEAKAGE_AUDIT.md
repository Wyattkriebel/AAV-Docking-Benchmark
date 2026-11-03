# Executed source-only family and training-overlap audit

**Seven of the 12 new structural candidates share an annotated CATH homologous
superfamily with a protein in the default P2Rank training dataset. Two share a
CATH homologous superfamily with development targets.** Zero exact PDB-ID
matches would have detected either problem. The audit certifies no universally
novel or learned-method-training-held-out families.

The [machine-readable audit](../benchmarks/family_overlap_audit.json) records
actual alignments, source annotations, warnings and missing evidence. All
[12 candidate identities](../benchmarks/heldout_candidate_inventory.json) remain
unchanged. No new docking or pocket-prediction outcomes were read. Structural
reference qualification and this family audit are separate gates.

## Sources, denominators and rule frozen before calculation

The default P2Rank model's pinned provenance identifies **CHEN11-fpocket**. Its
[official dataset list](https://raw.githubusercontent.com/rdk/p2rank-datasets/0236ecb38cbb60b89849a1ea36fbd9ee93f3e906/chen11-fpocket.ds)
has **251 rows and 241 unique PDB entries**, not 11 proteins. Downloaded bytes
match SHA256 `57344f47460fe786a62f16eec1168018b5d277008ced13dde206ecf4b2fd8419`
in [p2rank.lock.json](../benchmarks/p2rank.lock.json). There are no exact
candidate/training or development/training PDB-ID overlaps.

The source collector retrieved all **262 distinct PDB entries** across those
241 training entries, the 12 fixed candidates and nine development entries
(154L, 1GNY, 1KJL, 1QFO, 2J1T, 1PWB, 2VUZ, 7MTZ, 7MT0). Seventeen archived
[RCSB GraphQL](https://data.rcsb.org/) requests and responses include all
**297 deposited protein entities**, canonical sequences, author/label chain
maps, CATH/SCOP annotations and domain boundaries. There were no missing entries
or GraphQL errors. Current source metadata is pinned; it is not a claim that
these are precisely the historical training-coordinate revisions.

The [qualification policy](../data/archive/family_audit/qualification_policy.json)
was written and hashed before sequence similarities were calculated. It uses:

- Local affine-gap alignment with BioPython 1.86, BLOSUM62, gap opening −10,
  extension −0.5 and epsilon 10⁻⁶. The implementation reports the Gotoh local
  alignment algorithm. The first highest-score alignment is selected with
  deterministic record-ID ordering and all breakpoints retained.
- Identity ≥30% of aligned **nongap residue pairs**, and ≥70% coverage of the
  shorter compared sequence. Coverage counts residue pairs, excluding residues
  aligned to gaps. Unknown/ambiguous symbols cannot supply identities. These are
  an operational warning screen, not a significance-calibrated homology test.
- Exact four-level CATH **homologous-superfamily** matches as separate evidence.
  Shared three-level topology is separately recorded. Missing classification
  stays **unknown**. A novel PDB ID is never sufficient.

CATH groups individual domains; its third level describes topology/fold and its
fourth level concerns evolutionary homology. These are different levels of
independence. See the [CATH hierarchy](https://cathdb.info/wiki/doku/?id=faq)
and [homologous-superfamily definition](https://www.cathdb.info/wiki/doku/?id=glossary%3Ahomologous_superfamily).
The exact gap and local-alignment semantics follow the
[BioPython pairwise-alignment documentation](https://biopython.org/docs/latest/Tutorial/chapter_pairwise.html),
with the installed version and every scoring parameter explicitly pinned.

Whole-chain sequences and source-defined CATH/SCOP2 superfamily domains were
screened. Discontinuous segments with the same CATH domain ID are joined in
sequence order. Identical domain definitions across crystal copies are
deduplicated; differing source boundaries remain visible. This produced **1,365
sequence records**, including 42 candidate chain/domain records, 1,293 training
records and 30 development records. The complete **56,366 pairwise alignments**
cover candidates against development, training and other candidates. No
sequence-negative pair is omitted from the archived denominator.

## What the comparisons found

The table separates an exact CATH homologous-superfamily match from a numerical
sequence warning. A warning can involve a short segment, non-recognition domain
or gap-rich cross-fold alignment and therefore does not automatically establish
common ancestry. It also prevents an unqualified novelty assertion until
source-based adjudication resolves it. All individual hits and breakpoints are
available in `all_alignments.json.gz` and indexed from the audit JSON.

| Candidate | Development comparison | P2Rank training CATH match | Interpretation |
|---|---|---|---|
| 1I3H, ConA | CATH `2.60.120.200` matches galectin 1KJL; numerical warnings also present | 1URX, `2.60.120.200` | Shared development/training superfamily; new name does not supply independence. |
| 2G7C, TcdA CROPS | No CATH match; sequence warnings against AAV9 fragments | 2J8F, `2.10.270.10` | Training includes a matching repeat-domain superfamily; sequence independence from development remains unadjudicated. |
| 2GUP, ROK | No CATH match; sequence warnings against 2VUZ and a short 1PWB coiled-coil segment | 1CZA and 3BEX, `3.30.420.40` | Training superfamily seen. Numerical warnings are not a proof of lectin/kinase homology. |
| 2RDK, cyanovirin | No sequence warning or CATH match to development | 1IIY, `2.30.60.10` | Development overlap not detected; training contains cyanovirin itself. |
| 2WAB, CE2 | CATH `2.60.120.260` matches 1GNY/2J1T | 1UYY, `2.60.120.260` | Shared N-terminal domain blocks whole-protein family independence; recognition-domain context must be explicit. |
| 2WNB, GT29 | No sequence warning or CATH match to development | None detected | Training sequence warnings and classification gaps prevent an unseen-family claim. Full aglycone chemistry remains a separate gate. |
| 3A2S, PorB | No CATH match; warnings involve a short 1PWB coiled-coil domain | None detected | Warning credibility and missing comparison classifications remain unresolved. |
| 3ZWE, BambL | No sequence warning or CATH match to development | None detected | Training sequence warnings remain; exact PDB is also in ProtGlycanDock train. |
| 4R9F, CpMnBP1 | No sequence warning or CATH match to development | 1LBC, `3.40.190.10` | Development overlap not detected; training shares SBP-II superfamily. |
| 5HZA, norovirus P domain | No sequence warning or CATH match to development | None detected | Training sequence warnings/classification gaps remain; no certification of novelty. |
| 5V6F, RbmC β-prism | No CATH match; numerical warnings against AAV9 and a short 1PWB segment | 2HYR, `2.100.10.30` | Training contains griffithsin in the same annotated superfamily. |
| 8I5U, GH116 | CATH missing; sequence warning against 2J1T | CATH missing | SCOP2 six-hairpin/galactose-mutarotase annotations and training matches require domain-level interpretation; remains unknown, not novel. |

Three concrete examples show why both sequence and structural evidence matter:

- **2RDK versus training 1IIY:** 96 identities among 101 aligned pairs
  (**95.05% identity, 100% shorter-sequence coverage**), plus identical CATH
  `2.30.60.10`. The training source describes cyanovirin-N, and its author chain
  matches the chain hint in the official dataset path.
- **4R9F versus training 1LBC:** 61/204 identities (**29.90%**) narrowly misses
  the sequence cutoff, while both have CATH `3.40.190.10`. The latter is a
  glutamine-binding receptor; ligand function need not match for homologous
  superfamily exposure.
- **5V6F versus training 2HYR:** 35/89 identities (**39.33% identity, 72.95%
  shorter coverage**), plus identical CATH `2.100.10.30`. The training source is
  griffithsin. Its dataset-selected author chain is confirmed in the metadata.

The comparisons found sequence-screen warnings against development for seven
candidates; only **2RDK, 2WNB, 3ZWE, 4R9F and 5HZA** had neither a sequence
warning nor a CATH homologous-superfamily match in the available development
classifications. This is a finite screened-source result, not certification of
five independent novel folds.

## Unresolved evidence and permitted claims

CATH is present for 11/12 candidate entities, 7/9 development entities and
266/276 training entities. Development AAV9 holo/apo entries 7MTZ/7MT0 and
candidate 8I5U lack a CATH assignment in the pinned snapshots. Two
**dataset-selected** training entities, 1ZB6_1 and 2HXF_2, lack CATH. The other
unclassified training co-entities include short peptides. Their exact sequences
were still included in the alignment screen. Most CATH annotations in these
RCSB snapshots are from version 4.3; newer unclassified structures do not become
novel by absence.

The collector conservatively included every protein entity in a training PDB,
even if the official dataset path selects another chain. Each hit records
whether its entity matches that path's author-chain hint. All seven positive
CATH-overlap examples in the table involve an entity matching a selected hint;
they are not caused solely by an unrelated co-entity. Some numerical sequence
warnings do involve short co-entities or small domains and need interpretation.

Percent identity was frozen with a nongap-pair denominator. Long, gap-rich local
alignments can satisfy that screen across proteins with different CATH folds.
The audit therefore retains warnings without equating every threshold crossing
to credible homology. Gap-inclusive identity is supplied as a descriptive
quantity, with **no change to the frozen primary screen**. Resolving warning
credibility requires independent recognition-domain/structural evidence, not
inspection of docking success or changing cutoffs to increase a desired cohort.

The original six priority map-review cases remain useful candidates, but they
are not six cleared P2Rank-held-out families. Cyanovirin, SBP-II and the RbmC
β-prism have direct training-superfamily exposure; GH116 lacks primary CATH
classification and has SCOP2 superfamily matches to training entries 1V7W,
1W3Y and 2B4F. BambL and norovirus have numerical warnings and incomplete
comparison classification. A prospective external workflow check can still be
performed with the uncertainty and exposure strata stated accurately.

This audit covers the **declared default-model CHEN11-fpocket protein list**.
Feature-table construction, probability/score calibration, GNINA and other
learned methods remain separate provenance questions. No universal
learned-method-held-out claim is earned. No candidate inventory, docking input,
eligibility denominator or method parameter was changed based on outcomes.

## Archive and replay

The [artifact manifest](../data/archive/family_audit/artifact_manifest.json)
pins raw official dataset bytes, all GraphQL requests/responses, canonical
protein/domain sequence records, the frozen policy, software parameters and all
alignment breakpoints/metrics. The verifier in
[family_audit.py](../aav_benchmark/family_audit.py) never opens docking, pocket or
pose-metric results.

Default replay needs only the repository Python environment. It checks source
hashes and complete denominators, rederives sequences/CATH assignments from raw
responses, checks every domain sequence, and recalculates all 56,366 identity,
coverage and threshold decisions from the retained breakpoints:

```sh
python -m aav_benchmark.family_audit
python -m pytest -q tests/test_family_audit.py
```

To independently recalculate every optimal local affine-gap alignment, use an
isolated environment with the pinned dependencies; the repository dependency
set was not changed:

```sh
python3.11 -m venv /tmp/aav-family-replay
/tmp/aav-family-replay/bin/pip install -r data/archive/family_audit/requirements.txt
/tmp/aav-family-replay/bin/python -m aav_benchmark.family_audit --recompute-alignments
```

The tests include analytic controls for exact identities, insertions, coverage,
ambiguous residues, inclusive thresholds, invalid coordinates and the complete
source denominator. SHA replay demonstrates byte and calculation consistency;
it does not independently prove manual biological homology interpretation.
