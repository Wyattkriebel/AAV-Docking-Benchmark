# Independent check of the site and Gal9C results

Review completed 1 November 2026. A separate calculation read the frozen source
ledger, preparation records, execution commands, raw Vina output and logs. It
did not import the benchmark's preparation, engine, metric or report functions.
The machine-readable record is
[`independent_review.json`](../results/gal9_specificity/independent_review.json).

The check passed 1,086 assertions. These are artifact and arithmetic checks, not
new experimental measurements or a second docking experiment.

## Gal9C: the reported arithmetic is reproduced

All 12 frozen compounds have exactly the three specified seeds, giving 36
successful runs and 20 returned modes per run. Every run's receptor, ligand,
configuration, preparation, reference and pose fingerprints match the local
files. The numerical command options match the frozen experiment, and MODEL 1
energies agree with the raw log tables. Aggregation uses the mean of the three
MODEL 1 energies for each compound.

The 108 bundled prepared-input archives also pass independent compressed and
decompressed checksum/size checks; their restored bytes equal the reviewed
preparation files.

| Independently calculated quantity | Result |
| --- | ---: |
| Correct cross-class energy comparisons, with half credit for ties | 26/36 |
| AUROC | 0.7222222222 |
| Threshold-group average precision | 0.7842592593 |
| Correctly ordered source-matched epimer pairs | 5/6 |
| Matched-pair concordance | 0.8333333333 |
| Spearman correlation of mean energy with finite measured KD, six positives | 0.7142857143 |
| Heavy-atom-count AUROC / AP, either score orientation | 0.5 / 0.5 |
| Formal-charge AUROC / AP | 0.5 / 0.5 |

Average precision was calculated by grouping equal scores, advancing recall at
each distinct score, and weighting the precision after that group by its recall
increment. AUROC was calculated directly from all 36 measured/NB comparisons.
Spearman correlation used average ranks; simultaneous reversal to −energy and
pKD leaves its sign unchanged.

All 12 prepared SDF graphs match their frozen isomeric SMILES. Reassigning stereo
from SDF coordinates and, separately, from each of the 36 MODEL 1 coordinate
sets preserves those same graphs. Every top pose contains the exact expected
heavy-atom identities. All six source-matched pairs have the same connectivity,
formula, heavy-atom count and zero formal charge; their only stereochemical
difference is the arylamino-bearing sugar C3 center. This check does not establish
an experimentally correct docked pose.

The 8/21 pair is ordered incorrectly by mean energy. Pairs 4/17 and 6/19 have
small favorable energy differences, about 0.0893 and 0.0427 kcal/mol. These
differences are reported as descriptive score differences, without a claim of
experimental precision or a fitted classification threshold.

All six NB rows retain null numerical KD, null numerical lower bound and null
KD relation. Their 2 mM value remains the source's highest tested concentration.
The source ledger and preparation labels agree exactly. The journal and
supporting-information bytes are not redistributed, so this check verifies
their frozen curated transcription and declared source checksums, rather than
claiming a fresh local verification of all publication bytes.

## Site localization: denominators and summaries agree

The endpoint table contains all 27 method/condition combinations: eight bound
complexes and the additional paired AAV9 apo condition, each evaluated with
fpocket, P2Rank and the random surface baseline. Independently averaging the
eight bound conditions gives the following primary 4 Å endpoint:

| Method | Top 1 | Top 5 |
| --- | ---: | ---: |
| fpocket | 3/8 = 0.375 | 5/8 = 0.625 |
| P2Rank | 3/8 = 0.375 | 4/8 = 0.500 |
| Random surface, mean conditional null success | 0.04125 | 0.20750 |

Both corrected full-capsid fpocket runs have return code −5, no predictions and
zero success at all three distance thresholds. They remain in their eligible
denominators. All nine original fpocket reader failures are retained separately.
All nine P2Rank records completed, including the two conditions with zero
predictions. Each baseline condition retains exactly the 100 frozen seeds and
five patches per seed. These draws characterize a conditional computational
null; they are not 100 independent biological cases. The apo condition is a
paired condition of AAV9 and is not added to the eight-complex primary mean.

The final evaluation refresh preserves all 27 primary 4 Å endpoint pairs and
all primary bound-condition macro means. Both prescribed 2 Å and 6 Å candidate
merge sensitivities are complete for every method/condition, giving 54
additional summaries. Their trial counts, endpoint aggregation and ordered
3/4/5 Å success thresholds were checked independently from the recorded
evaluation details. The JSON review records the final evaluation checksum.

This independent check aggregates the recorded site endpoints and checks run
statuses and baseline budgets. Raw pocket parsing, archived receptor restoration
and geometric endpoint re-evaluation are separate verification steps, required
in addition to this arithmetic check.

## Interpretation and a frozen-configuration clarification

The Gal9C experiment is a known-site comparison in one assay context and one
related synthetic glycoside series. The source-derived construct and tag
provenance remains approximate. No independent receptor-family test, calibrated
affinity prediction or atomic pose reference is established by these results.
The site panel remains a set of qualified descriptive published-site references;
successful localization does not establish binding specificity or atomic pose
accuracy. Compatibility retries and the source-only symmetry amendment must
remain visible when describing the experiment's development history.

The immutable Gal9C configuration inherited dormant pose thresholds and the
old `specificity` / `confirmatory_status` strings from its parent configuration.
The companion
[`gal9_field_interpretation.json`](../benchmarks/gal9_field_interpretation.json)
now identifies those fields as inactive and binds its explanation to the
unchanged frozen configuration hash, resolving this interpretation issue.
Active Gal9C settings are the frozen cohort, known site, seeds, docking budget,
energy aggregation and named specificity metrics. Input conformers are
explicitly ineligible pose references throughout the prepared records.
