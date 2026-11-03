# Independent review of the POL increment

Review date: 3 November 2026. A separate agent reviewed the complete frozen
experiment and archives without changing the preparation, raw evidence or
analysis. This is an implementation and interpretation review, not independent
biological replication or external preregistration.

- All 43 targeted preparation, raw-report and execution-label tests passed.
- The report verifier passed. A fresh archive-only relocation, without
  `data/processed/` or installed docking tools, reproduced exact JSON and Markdown.
- Independent parsing of all 96 raw Top-1 score records reproduced the primary
  galactose mean -4.2045, glucose mean -4.1771666667 kcal/mol and six of seven
  sensitivity reversals.
- All seven frozen preparation-helper hashes, nine execution snapshots and
  460 evidence artifacts verified.
- A fresh retrieval of the primary article matched the recorded SHA256 and
  byte count. Its assay text confirms actually tested glucose with no detected
  ITC binding, the approximate galactose affinity and the recorded buffer,
  injection and concentration context. Article bytes are not redistributed.
- Two experimental reducing sugars remain distinct from four modeled anomers.
  Equal anomer weighting, known-site information and categorical no-binding
  labels remain explicit. No classifier, affinity calibration or population
  reliability claim is made.

No material publishing blocker was identified. See the
[results](POL_SPECIFICITY_RESULTS.md), [source qualification](INDEPENDENT_SPECIFICITY_CANDIDATES.md)
and [raw-data summary](../results/pol_specificity/summary.json).
