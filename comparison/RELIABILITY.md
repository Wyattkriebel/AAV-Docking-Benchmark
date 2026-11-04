# What the evidence says about glycan-docking reliability

Assessment date: 3 November 2026.

**Docking has conditional predictive value, but a favorable docking score does
not establish biologically correct recognition. General reliability is not
established by this repository.** Its released nonviral development results
show partial site-localization utility, failed stringent pose recovery in the
tested workflows, and limited or preparation-dependent specificity ranking.
Published studies establish useful performance under particular conditions;
they do not turn these results into a universal success or failure claim.

This is a completed assessment of existing evidence. It is **not a completed
confirmatory benchmark**, new experimental validation, or a certification that
a method meets a useful-performance threshold. Those distinctions are recorded
in the [machine-readable evidence](reliability_evidence.json).

## Three questions, three different endpoints

| Question | Relevant evidence | What the available results support |
| --- | --- | --- |
| Does a workflow locate the modeled site? | Ranked pocket centers compared with a structural site reference | Partial recovery in the released nonviral cases; localization does not establish ligand pose or binding specificity. |
| Does it return and rank the modeled pose? | Chemically mapped ligand coordinates, evaluated in the receptor frame | The tested Vina workflows failed their primary 2 Å endpoint on the two full-model and four core-reference cases. |
| Does its score distinguish measured binding alternatives? | Tested alternatives within the same target and assay context | Imperfect ranking in one Gal9C series and an unstable score preference in one POL contrast; no general classifier or affinity calibration. |

These endpoints cannot be pooled into one reliability percentage. fpocket and
P2Rank are pocket predictors here; the pose experiments use Vina downstream.
They are not three independent docking algorithms.

A reliability claim needs a named workflow, endpoint and target population,
with uncertainty based on independent biological units. Practical use can also
require an agreed minimum performance. These development results do not test
such a population-level floor; the conclusion "unestablished" is not a
statistical demonstration that every method performs below a useful level.

## Repository evidence

The [evidence packet](reliability_evidence.json) contains selected fields from
four previously released, checksum-identified result files. Its source URLs
point to the retained research snapshot, not to an unreported new experiment.

### Site localization: a retrospective nonviral subset

Seven cases are included: 154L, 1GNY, 1KJL, 1QFO, 2J1T, 1PWB and 2VUZ.
The criterion is a candidate center within 4 Å of a modeled ligand heavy atom.
Top5 means five distinct candidate-site classes, not five docked ligand poses.

| Case | fpocket Top5 | P2Rank Top5 |
| --- | ---: | ---: |
| 154L | yes | yes |
| 1GNY | yes | yes |
| 1KJL | yes | yes |
| 1QFO | no | no |
| 2J1T | yes | no |
| 1PWB | no | no |
| 2VUZ | yes | yes |
| Descriptive total | 5/7 | 4/7 |

This nonviral slice was selected retrospectively from existing results. The
original eight-system primary result was 5/8 and 4/8, respectively; those
denominators are preserved in the
[original site report](https://github.com/AlecKriebel/AAV-Docking-Benchmark/blob/ddd84f9c9021d43a4c58868c9116d0523f3cbb79/docs/SITE_LOCALIZATION_RESULTS.md).
The seven-case fractions are not new prospective endpoints or estimates of
population reliability. References have local-density, chemistry and assembly
cautions; localization near a modeled ligand is not experimental confirmation
that a predicted pocket is a biologically functional site. Training overlap
and family independence have not been established for this subset.

### Pose recovery: distinguish full models from partial cores

The corrected V2 evaluation uses two complete modeled ligands (1GNY and 1PWB)
and four modeled cores (1KJL, 1QFO, 2J1T and 2VUZ). The strata remain separate.

| Vina input information | Full-model Top1 / Top5 / best retained at 2 Å | Core Top1 / Top5 / best retained at 2 Å |
| --- | --- | --- |
| fpocket-guided, five-site budget | 0% / 0% / 0% | 0% / 0% / 0% |
| P2Rank-guided, five-site budget | 0% / 0% / 0% | 0% / 0% / 0% |
| Known-site oracle, one box | 0% / 0% / 0% | 0% / 0% / 0% |

Each percentage averages three seeds within each complex, then weights complexes
equally within its stratum. Global and local-interface receptor alignments give
the same zero primary recovery; the one-site predicted-budget sensitivity is
also zero. Unsupported 1QFO preparation remains an eligible workflow zero.
No core result is relabeled as full-ligand recovery.

At the predefined 3 Å sensitivity, the 11-atom 1KJL core is recovered among Top5
and best-retained poses, but never Top1. This gives a core-stratum mean of 25%
for those endpoints, with no full-ligand success. The threshold sensitivity
matters and is reported rather than suppressed.

These are failures of the tested preparation-and-docking workflows under their
recorded limits. They do not isolate an intrinsic algorithm defect. Full ring
pucker and torsion-energy validity were not certified, and an absent near-native
returned pose does not reveal whether one was generated internally. V2 also
corrects an atom-correspondence error after diagnostic outcomes had been seen;
it is development evidence rather than pristine prospective confirmation.
See the [corrected source report](https://github.com/AlecKriebel/AAV-Docking-Benchmark/blob/ddd84f9c9021d43a4c58868c9116d0523f3cbb79/docs/QUALIFIED_BLIND_POSE_V2_RESULTS.md).

### Specificity: two assay contexts, different limitations

| Evidence | Observed result | Limit on the conclusion |
| --- | --- | --- |
| Gal9C: six measured binders and six tested matched epimers | AUROC 0.7222; average precision 0.7843; measured binder favored in 5/6 pairs. Size and charge baselines each give AUROC 0.5. | One target and synthetic series; no calibrated decision threshold, cross-family validation or quantitative affinity claim. |
| POL/CBM67: one measured galactose/glucose contrast | Primary glucose-minus-galactose score difference +0.0273 kcal/mol; direction reverses in 6/7 prespecified sensitivities. | One experimental contrast with computational repeats; the favorable primary sign is not a robust specificity result. |

Gal9C's one discordant matched pair supplies a counterexample to treating lower
Vina energy as a guarantee of the measured binding preference in that assay.
It does not establish a general error rate. POL further demonstrates that a
small score preference can depend on preparation and aggregation choices.

No-binding annotations are assay-conditional categorical observations; they
are not numerical dissociation-constant bounds. Docking seeds, modeled anomers,
site choices and preparation sensitivities are not independent experimental
replications. Exact size/charge matching removes those two simple separators
within the selected Gal9C panel; it does not demonstrate that matching itself
improved predictive performance. The original
[Gal9C](https://github.com/AlecKriebel/AAV-Docking-Benchmark/blob/ddd84f9c9021d43a4c58868c9116d0523f3cbb79/docs/GAL9_SPECIFICITY_RESULTS.md)
and [POL](https://github.com/AlecKriebel/AAV-Docking-Benchmark/blob/ddd84f9c9021d43a4c58868c9116d0523f3cbb79/docs/POL_SPECIFICITY_RESULTS.md)
reports retain the experimental provenance and model limitations.

## How this fits the primary literature

This is a targeted comparison with relevant primary studies, not a systematic
review or a meta-analysis. These reported results use different inputs,
budgets, endpoints and test sets; their percentages must not be ranked or pooled
as if they came from one controlled experiment.

| Primary study | Reported evidence | Interpretation boundary |
| --- | --- | --- |
| [Ranaudo et al., HADDOCK, 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11480977/) | 89 bound complexes and a 55-case unbound subset; approximately 70% and 40% Top5 success reported in the abstract with binding-site information. | Supports conditional pose-modeling utility. The study uses interface-ligand RMSD and several quality cutoffs; these headline rates are not a blind-site or binder/nonbinder result. |
| [Nance et al., GlycanDock, 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8742512/) | Local refinement evaluated on 109 bound and 62 unbound structures, with controlled starting deviations and a ring-atom success statistic. | Refinement from a putative complex is a different task from locating an unknown site. The paper's one specificity example is not a general validated classifier. |
| [Xu et al., ProtGlycanDock, October 2026 preprint](https://arxiv.org/html/2610.26009v1) | Ten methods evaluated with protein- and glycan-holdout tests of 148 and 134 samples. | Adds substantially broader structural evidence. Structures with more than four protein chains are excluded, and its docking baselines use predicted receptors and pockets. These conditions differ from this repository's oracle and experimental-receptor cases. |
| [Zhao et al., GlycoWeave, October 2026 preprint](https://arxiv.org/abs/2610.37410v1) | The abstract reports 71.9% Top25 high-quality success on 698 complexes. | An author-reported preprint result with a Top25 endpoint; it is not a Top1 guarantee, assay-specificity estimate or replication by this repository. |

The inference from these studies is conditional utility, not universal
incapability or universal reliability. Structural prediction success also does
not by itself establish experimental binding discrimination.

## What remains unresolved

The evidence does not establish a useful-performance floor on an independently
held-out population of protein families. The small development collection,
related folds, unresolved training overlap and preparation changes prevent
that inference. No family-level reliability interval is reported here, and
seeds or repeated poses cannot supply the missing biological replication.

There is also no independently validated score threshold for sensitivity and
specificity, no calibrated affinity predictor, and no common new-cohort
comparison of multiple docking algorithms. Reproducing a released example
validates execution; it does not fill those gaps. The wider confirmatory study
in [PROTOCOL.md](../PROTOCOL.md) therefore remains incomplete.

## Audit and reproduction boundaries

The evidence packet identifies source paths, immutable commit URLs and SHA256
digests. Pose aggregates are copied from corrected V2 without changing its
strata or estimands; site rows are explicitly labeled as a retrospective slice.
Gal9C and POL remain separate assay contexts. Exact copied-field equality and
the stated descriptive arithmetic were checked against the retained sources.

The existing `python3 pilot/pilot-script.py --check` command reproduces the
separate six-run AAV9 comparison. It does **not** independently reproduce this
broader evidence assessment or the cited studies. This assessment adds no
native docking, atomic scoring, a newly curated cohort or new biological assay
data. The later mixed cohort and supplemental expansion are not evaluated or
silently filtered into a new benchmark.
