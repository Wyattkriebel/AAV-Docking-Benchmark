# Independent source-cohort PR review

Review date: 3 November 2026. This review inspected source qualification,
assembly-copy truth, sequence/domain overlap, archive portability and publication
packaging. It did not inspect new test pocket predictions, docking poses, scores
or recovery endpoints.

The numerical evidence is ready for a source-only PR. Publication packaging
requires the eleven exclusions enumerated in
[`source_public_distribution.json`](../benchmarks/source_public_distribution.json).
Those local files remain unchanged; their hashes, sizes, retrieval URLs and
reviewed access scope remain in the original provenance. They must not be added
to the public Git tree. The separate public verifier permits exactly those
noncomputational snapshots to be absent and keeps every numerical source check.

## Reference and assembly findings

The fixed twelve-case ledger reproduces six full modeled ligands, five whole
residue cores and one atomic rejection. It retains 388 source atom rows,
377 unique parent heavy atoms, 298 primary scored atoms, 19 same-chemical
asymmetric-unit copies and 16 qualified biological-assembly instances. The
six-case starting set is fixed as 2RDK, 3ZWE, 4R9F, 5HZA, 5V6F and 8I5U:
five full modeled ligands and the 11-atom BMA B2 core of 4R9F. Eight candidates
meet source-conditioned preparation requirements; that count does not establish
method compatibility or execution success.

The replay independently reconstructs actual deposited coefficient maps,
residue gates, canonical masks and biological-assembly transforms. The maximum
observed map-sample difference on this host is zero sigma. Copies must pass the
same primary canonical mask and belong to the selected assembly. The valid
2RDK ligand in another assembly cannot be substituted; the failing full 3ZWE
copy and its chemically different partial sites remain explicit unscored known
sites. A core never authorizes truncating the preparation parent. In particular,
2WNB retains the full 35-atom sugar/aglycone graph despite its 25-atom scoring
core, and 2WAB retains its unresolved parent alternatives as a preparation block.

No source replay proves a unique ligand conformation, removes model-phased map
bias or independently re-establishes the curator's paper/reagent/function
judgments. Residue-average `avgoccu` is not an all-atom occupancy minimum. These
limits are correctly disclosed in the qualification report.

## Family and training findings

The source audit covers the official default-model CHEN11-fpocket list: 251 rows,
241 unique PDB entries. All 262 combined entries, 297 protein entities,
1,365 chain/domain sequence records and 56,366 candidate comparisons reproduce
from pinned source responses and alignment breakpoints. Two candidates share a
CATH homologous superfamily with development targets; seven share one with a
dataset-selected P2Rank training entity. Zero exact PDB-ID overlaps would miss
those relationships. Missing classifications and numerical warnings remain
unknown rather than certified novelty. The five screen-negative development
comparisons are not five proven independent folds.

The earlier inventory's suggested 30% identity/80% coverage screen was a proposal.
The actual separately frozen audit policy uses 30% nongap-pair identity and 70%
shorter-sequence coverage, before any similarity calculation. That distinction
must remain visible. Default replay recomputes all archived identity, coverage
and threshold decisions; optimal affine-gap alignment recalculation is a
separate BioPython 1.86 option. Protein-list exposure does not resolve feature
construction, probability calibration, GNINA or other model provenance.

## Publication redistribution correction

The [PMC copyright notice](https://pmc.ncbi.nlm.nih.gov/about/copyright/)
distinguishes free access from permission to redistribute, and directs users to
the individual article's terms. PDB archive CC0 does not cover publisher text or
figures. The original source manifest's statement that papers have separate
copyright contexts is accurate but insufficient to authorize redistribution.

The public distribution omits these exact original snapshots:

| Snapshot basename under `data/archive/heldout_reference/papers/` | Reason |
|---|---|
| `1I3H.primary_abstract.html` | Publisher-copyright abstract; no redistribution license established. |
| `2G7C.primary_abstract.html` | Publisher-copyright abstract; no redistribution license established. |
| `2GUP.deposition_page.html` | Provider web-page snapshot; numerical PDB/API records remain available. |
| `2RDK.PMC2327273.html` | Protein Society copyright; archived Crossref points to Wiley terms rather than an open reuse license. |
| `2WNB.publisher_response.html` | Closed publisher page/response; no open redistribution license established. |
| `2WNB.supplement.pdf` | Publisher supplement to that closed article; no open redistribution license established. |
| `3A2S.PMC2872391.html` | No explicit open redistribution license in the archived article. |
| `3ZWE.supplement.pdf.response.html` | Browser-check response, not a recovered scientific supplement. |
| `5HZA.PMC4836343.html` | ASM copyright notice states all rights reserved. |
| `5HZA.primary_density_figure.jpg` | Figure from the same uncleared article. |
| `8I5U.supplement.pdf.response.html` | Browser-check response, not a recovered scientific supplement. |

The original frozen manifest remains byte-exact with 174 records. The public
tree contains 163 of those archives, with eleven metadata-only omissions. No
coordinate, validation, CCD, structure factor, map coefficient, scored mask,
assembly operator, family source or alignment is omitted. All available optional
local snapshots still must match their original hashes. Public replay does not
download or create placeholder articles and does not claim to re-prove manual
publication judgments.

Retained papers and figures have separate attribution, source and license records
in the new distribution manifest. Their bytes remain unmodified and do not take
the repository's software license:

| Source | Retained article/figure | License basis |
|---|---|---|
| Montanier et al., 2009, [doi:10.1371/journal.pbio.1000071](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1000071) | 2WAB article | Explicit Creative Commons Attribution notice in both archived article and publisher page. The notice does not state a version; none is invented. |
| Audfray et al., 2012, [doi:10.1074/jbc.M111.314831](https://doi.org/10.1074/jbc.M111.314831) | 3ZWE article and Figure 4 | Publisher-deposited [Crossref metadata](https://api.crossref.org/works/10.1074/jbc.M111.314831) specifies CC-BY-4.0 for the version of record from 10 September 2020; archived metadata and fresh retrieval agree. |
| Chekan et al., 2014, [doi:10.1074/jbc.M114.579904](https://doi.org/10.1074/jbc.M114.579904) | 4R9F article and Figure 4 | Publisher-deposited [Crossref metadata](https://api.crossref.org/works/10.1074/jbc.M114.579904) specifies CC-BY-4.0 for the version of record from 10 October 2020; archived metadata and fresh retrieval agree. |
| De et al., 2018, [doi:10.1371/journal.ppat.1006841](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1006841) | 5V6F article and Figure 4 | Explicit CC-BY-4.0 notice in archived article and matching publisher notice. |
| Pengthaisong et al., 2023, [doi:10.1021/acscatal.3c00620](https://doi.org/10.1021/acscatal.3c00620) | 8I5U article | Archived ACS article explicitly links CC-BY-4.0 and attributes the authors. |

## Documentation corrections and validation

`HELDOUT_COHORT_DESIGN.md` now identifies its zero-acceptance inventory fields and
availability table as the earlier snapshot, links the completed qualification
and overlap audit, distinguishes the seven map-review priorities from the fixed
six-case starting set, and records the proposal-versus-executed screen thresholds.
The fixed JSON inventory remains unchanged. The qualification report already
corrects the inventory's erroneous *Clostridium perfringens* expansion of 4R9F:
the deposited source and primary paper identify *Caldanaerobius
polysaccharolyticus*. Readers must use that explicit erratum.

The qualification report's original phrases “full text was archived” and “figures
were archived” describe the local research archive. The source PR should link
this public distribution manifest there so those phrases cannot imply that every
retrieved article is redistributed. Its public reproduction command must use:

```sh
python scripts/verify_public_source_cohort.py
python -m pytest -q tests/test_public_source_replay.py tests/test_heldout_reference.py tests/test_family_audit.py
```

All **33 targeted tests pass**, including a real numerical replay with all eleven
snapshots absent. Controls reject added omissions, map loss, changed original
manifests, rebound omission hashes, corrupt retained or optional article bytes,
and attempts to consume waived publication bytes numerically. The original
strict verifier still rejects absent publication bytes. Numerical mutation
controls use private synthetic publication-free fixture manifests; no frozen
scientific manifest is rewritten. The hook is restored after failures.

Both source and family verifiers also pass from an archive-only relocated tree
with no results, processed inputs or executable tools. The new combined public
CLI additionally passes from a fresh relocated tree with all eleven excluded
snapshots absent: 163 retained records, 16 qualified assembly copies and 56,366
family comparisons verify. These validate source portability and arithmetic,
not docking performance, assay discrimination or completed scientific closure.
