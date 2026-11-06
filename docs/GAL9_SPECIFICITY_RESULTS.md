# Gal9C matched-epimer specificity results

This fixed analysis compares six compounds with finite measured affinities and six matched epimers reported as NB at the highest tested concentration, 2 mM. NB remains a categorical assay outcome; it is not a numerical KD bound.

The fixed denominator is 12 compounds × 3 seeds = 36 runs, covering six matched pairs. Verified successful runs: 36/36; failed: 0; running: 0; missing: 0. Complete compounds: 12/12.

AutoDock Vina 1.2.7 uses CPU 2, exhaustiveness 16, up to 20 modes and the same 24 Å box for every compound. Seeds are 104729, 130363 and 155921. Each compound's score is the mean of its three seed-specific top-ranked energies. Lower energy predicts measured binding. The site comes from the independent prior 3NV3 glycan structure; this is a known-site experiment.

Descriptive AUROC is 0.7222; step-wise average precision is 0.7843 at measured-binder prevalence 0.5000. No score threshold or classifier was fitted.

Across six source-matched epimer pairs, measured binders have lower mean energy in 5, conditional NB epimers have lower mean energy in 1, and 0 are exact aggregate ties. Ties receive half credit; concordance is 0.8333.

Among the six finite-affinity positives only, Spearman correlation of −mean Vina energy with measured pKD is 0.7143. NB compounds contribute no quantitative affinity value, correlation or censored numeric contrast.

| Descriptor baseline | AUROC | Average precision |
| --- | --- | --- |
| heavy_atom_count | 0.5000 | 0.5000 |
| formal_charge | 0.5000 | 0.5000 |

Heavy-atom-count AUROC is also 0.5000 in the reverse orientation because each matched pair has equal size. Formal charges are all zero.

| Compound | Series | Aryl group | Source assay result | Seed top-1 energies (kcal/mol) | Mean energy (kcal/mol) | Successful seeds |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | galacto | phenyl | 678 ± 30 µM | -6.206, -6.215, -6.225 | -6.2153 | 3/3 |
| 14 | gulo | phenyl | NB at 2 mM | -5.545, -5.530, -5.503 | -5.5260 | 3/3 |
| 4 | galacto | 4-methylphenyl | 770 ± 25 µM | -5.537, -5.542, -5.593 | -5.5573 | 3/3 |
| 17 | gulo | 4-methylphenyl | NB at 2 mM | -5.477, -5.467, -5.460 | -5.4680 | 3/3 |
| 5 | galacto | 3-methylphenyl | 670 ± 23 µM | -6.623, -6.625, -6.628 | -6.6253 | 3/3 |
| 18 | gulo | 3-methylphenyl | NB at 2 mM | -5.953, -5.953, -5.948 | -5.9513 | 3/3 |
| 6 | galacto | 3,4-dimethylphenyl | 130 ± 9 µM | -6.638, -6.649, -6.607 | -6.6313 | 3/3 |
| 19 | gulo | 3,4-dimethylphenyl | NB at 2 mM | -6.598, -6.569, -6.599 | -6.5887 | 3/3 |
| 7 | galacto | 3-chlorophenyl | 160 ± 6 µM | -6.380, -6.365, -6.363 | -6.3693 | 3/3 |
| 20 | gulo | 3-chlorophenyl | NB at 2 mM | -6.224, -6.201, -6.218 | -6.2143 | 3/3 |
| 8 | galacto | 3,4-difluorophenyl | 220 ± 7 µM | -5.652, -5.616, -5.579 | -5.6157 | 3/3 |
| 21 | gulo | 3,4-difluorophenyl | NB at 2 mM | -5.622, -5.705, -5.648 | -5.6583 | 3/3 |

Finite affinity errors above are source-reported SD from 4–8 single/double-point measurements. Seed energies follow the frozen seed order; incomplete compounds have no mean.

| Measured/NB epimer pair | Measured minus NB mean energy (kcal/mol) | Concordance |
| --- | --- | --- |
| 1/14 | -0.6893 | 1.0 |
| 4/17 | -0.0893 | 1.0 |
| 5/18 | -0.6740 | 1.0 |
| 6/19 | -0.0427 | 1.0 |
| 7/20 | -0.1550 | 1.0 |
| 8/21 | 0.0427 | 0.0 |

The 3NV3 declared 138-residue structural core exactly matches the corresponding assay sequence, including the annotated binding region. The model contains 136 observed residues: 13 N-terminal assay residues lie outside the declared structural core and two additional declared residues are unmodeled. Fusion/batch provenance remains approximate. One independently generated conformer is used per compound. No atomic pose truth for these alternatives is established and no input-coordinate RMSD is reported.

The frozen configuration retains unused fields inherited from the earlier structural diagnostic. The [field interpretation](../benchmarks/gal9_field_interpretation.json) identifies them without changing its frozen bytes. The explicit Gal9 panel, score aggregation and specificity metrics govern this analysis; pose thresholds are inactive.

These results describe one Gal9C assay context and one synthetic glycoside series. Gal9 and Gal3 share a fold; this provides no independent receptor-family generalization, calibrated affinity prediction or reliability certification.

All aggregates are recomputed from raw pose scores after checking the frozen experiment, source receptor, preparation/reference/ligand/receptor hashes, command options, raw pose hashes, full heavy-atom maps and raw log score tables. The JSON summary lists every expected run and a source hash manifest. Published assay-source checksums are preserved as declarations; publication bytes are not redistributed or claimed locally reverified.

Assay source: [Mahanti et al., Table 1 and supporting methods](https://doi.org/10.1021/acsmedchemlett.9b00396). Receptor qualification: [Gal9C receptor audit](../benchmarks/gal9_receptor_audit.json). Frozen protocol: [gal9_experiment.json](../benchmarks/gal9_experiment.json). Machine-readable evidence: [summary.json](../results/gal9_specificity/summary.json).
