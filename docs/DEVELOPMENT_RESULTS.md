# Development docking results

Generated deterministically by `python scripts/regenerate_report.py` from the archived summaries, run manifests, configuration files, per-pose metrics and pose artifacts listed below. `--check` verifies that this document is current. No raw result is modified.

**The repository's original scientific purpose remains incomplete.** These are AutoDock Vina 1.2.7 development diagnostics using true-site information. They do not establish blind site recovery, held-out pose reliability, discrimination of tested nonbinders, or reliability of other docking methods. No confirmatory cohort has been evaluated here.

## Eligible cases and preparation

The frozen external development selection retains **7 eligible cases**: 6 successful preparations, 1 failed and 0 unsupported; 0 still pending/running. These are preparation outcomes, not molecular binding outcomes. Docking accuracy below is conditional on available completed preparation and execution; absent docking is explicit.

| Case | Preparation | Docking evidence | Failure category | Reason |
| --- | --- | --- | --- | --- |
| 154L | ok | summary present | — | — |
| 1GNY | ok | summary present | — | — |
| 1KJL | ok | summary present | — | — |
| 1PWB | ok | summary present | — | — |
| 1QFO | failed | not completed / no summary | chemistry_reference_mismatch | CCD/reference stereochemistry mismatch 1:C1 |
| 2J1T | ok | summary present | — | — |
| 2VUZ | ok | summary present | — | — |

- **1PWB**: Native [Ca+2], fixed +2 Vina Ca atom; no specialized coordination/polarization.
- **2VUZ**: Resolved seven-residue glycan fragment; not whole nonasaccharide atomic truth; Native [Ca+2], fixed +2 Vina Ca atom; no specialized coordination/polarization; Missing depositedauthor91 restored to topology via audited numbering map; no loop atoms or cross-gap bond added; source outputidentities retained.

AAV9 is a separate diagnostic. Its apo/holo and complete-neighborhood/monomer conditions are paired model conditions for one biological system, not four independent cases. Monomer conditions are deliberately incomplete controls. External released receptors are not certified complete biological assemblies. Source reproduction alone does not qualify local density, stereochemistry, assemblies or family independence for confirmatory inference.

External protein-group labels are provisional descriptive annotations. In particular, 1PWB and 2VUZ share the broad C-type lectin fold; no count of independently qualified families or across-family generalization follows from this selection.

## Structural pose diagnostics

Each completed condition must contain its three configured seeds. Values are means across completed seeds within a case and condition; Top5 is the lowest RMSD among the five highest-ranked poses, and best retained is the lowest RMSD across all returned poses. Ligands are evaluated in the receptor coordinate frame without independent ligand alignment. No across-family confidence interval is estimated from these unqualified development cases.

| Case | Condition | Completed/eligible seeds | Mean Top1 RMSD Å | Mean Top5 minimum RMSD Å | Mean best retained RMSD Å |
| --- | --- | --- | --- | --- | --- |
| 154L | holo_complete | 3/3 | 8.890 | 7.478 | 2.446 |
| 1GNY | holo_complete | 3/3 | 14.249 | 6.458 | 3.661 |
| 1KJL | holo_complete | 3/3 | 7.071 | 6.488 | 4.861 |
| 1PWB | holo_complete | 3/3 | 8.024 | 4.432 | 3.861 |
| 2J1T | holo_complete | 3/3 | 7.665 | 5.580 | 5.238 |
| 2VUZ | holo_complete | 3/3 | 7.130 | 5.692 | 5.422 |
| AAV9 | apo_complete | 3/3 | 9.055 | 5.098 | 3.770 |
| AAV9 | apo_monomer | 3/3 | 9.683 | 8.999 | 6.440 |
| AAV9 | holo_complete | 3/3 | 11.494 | 8.470 | 3.333 |
| AAV9 | holo_monomer | 3/3 | 10.071 | 8.649 | 7.568 |

Counts below retain all configured seeds in each denominator. Repeated seeds and poses do not supply independent biological replication. A near-native pose absent from retained outputs does not identify internal sampling failure.

| Case | Condition | Top1 ≤2 Å | Top5 ≤2 Å | Retained ≤2 Å | Top1 ≤3 Å | Top5 ≤3 Å | Retained ≤3 Å |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 154L | holo_complete | 0/3 | 0/3 | 0/3 | 0/3 | 0/3 | 3/3 |
| 1GNY | holo_complete | 0/3 | 0/3 | 0/3 | 0/3 | 0/3 | 0/3 |
| 1KJL | holo_complete | 0/3 | 0/3 | 0/3 | 0/3 | 0/3 | 0/3 |
| 1PWB | holo_complete | 0/3 | 0/3 | 0/3 | 0/3 | 0/3 | 0/3 |
| 2J1T | holo_complete | 0/3 | 0/3 | 0/3 | 0/3 | 0/3 | 0/3 |
| 2VUZ | holo_complete | 0/3 | 0/3 | 0/3 | 0/3 | 0/3 | 0/3 |
| AAV9 | apo_complete | 0/3 | 0/3 | 0/3 | 0/3 | 0/3 | 0/3 |
| AAV9 | apo_monomer | 0/3 | 0/3 | 0/3 | 0/3 | 0/3 | 0/3 |
| AAV9 | holo_complete | 0/3 | 0/3 | 0/3 | 0/3 | 0/3 | 0/3 |
| AAV9 | holo_monomer | 0/3 | 0/3 | 0/3 | 0/3 | 0/3 | 0/3 |

## Contacts and physical overlap

Native heavy-atom contact identities use the fixed holo reference and the <4 Å cutoff. The generator checks that native contact and interface-group denominators remain constant across seeds and conditions for each case. Contact values are seed means for Top1. Interface groups follow receptor PDB residue records and include retained cofactors such as CA2; these values do not measure protein-residue-only coverage. Severe overlap means at least one evaluated ligand–receptor atom pair below 1.5 Å; this distance screen does not validate ring geometry, linkages or chemical plausibility.

| Case | Condition | Native atom contacts | Native interface groups | Contact recall | Contact precision | Interface-group recall | Interface-group precision | Top1 overlap seeds | Overlap evaluation scope |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 154L | holo_complete | 62 | 12 | 0.022 | 0.013 | 0.861 | 0.628 | 0/3 | provided_receptor_only |
| 1GNY | holo_complete | 56 | 9 | 0.036 | 0.023 | 0.889 | 0.533 | 0/3 | provided_receptor_only |
| 1KJL | holo_complete | 47 | 7 | 0.170 | 0.142 | 0.429 | 0.375 | 0/3 | provided_receptor_only |
| 1PWB | holo_complete | 37 | 7 | 0.108 | 0.066 | 0.714 | 0.556 | 0/3 | provided_receptor_only |
| 2J1T | holo_complete | 45 | 8 | 0.052 | 0.035 | 0.792 | 0.704 | 0/3 | provided_receptor_only |
| 2VUZ | holo_complete | 72 | 13 | 0.083 | 0.064 | 1.000 | 0.766 | 0/3 | provided_receptor_only |
| AAV9 | apo_complete | 35 | 3 | 0.029 | 0.021 | 0.333 | 0.167 | 0/3 | full_biological_assembly |
| AAV9 | apo_monomer | 35 | 3 | 0.019 | 0.021 | 0.333 | 0.178 | 3/3 | full_biological_assembly |
| AAV9 | holo_complete | 35 | 3 | 0.000 | 0.000 | 0.222 | 0.083 | 0/3 | full_biological_assembly |
| AAV9 | holo_monomer | 35 | 3 | 0.057 | 0.053 | 0.333 | 0.167 | 3/3 | full_biological_assembly |

## Native-coordinate scoring diagnostic

Experimental native heavy coordinates were supplied directly; RDKit added modeled hydrogens without optimization. These four score-only evaluations are diagnostics, **never docking successes**, and are excluded from the docking run/pose totals. Positive native-minus-docked gaps mean the prepared Vina score favors the returned Top1 pose over the supplied reference. This is a score-consistency observation for these inputs, not evidence about internal search sampling or measured affinity.

| AAV9 condition | Native heavy RMSD Å | Native score kcal/mol | Mean docked Top1 kcal/mol | Native minus docked kcal/mol | Score-only wall s |
| --- | --- | --- | --- | --- | --- |
| apo_complete | 0.000000 | 8.663 | -4.617 | 13.280 | 0.427 |
| apo_monomer | 0.000000 | 9.417 | -4.527 | 13.944 | 0.326 |
| holo_complete | 0.000000 | -1.858 | -4.866 | 3.008 | 0.455 |
| holo_monomer | 0.000000 | -1.037 | -4.120 | 3.083 | 0.328 |

## Galectin-3 quantitative assay panel

This exploratory panel contains three measured binders from one ITC context. All three compounds bind; there are no tested nonbinders, so binary AUROC/AUPRC, sensitivity and specificity are ineligible. KD values are published fitted estimates (±SEM), not raw ITC data; their uncertainty is not propagated into these descriptive rank correlations. The docking receptor is the 1KJL carbohydrate-recognition domain whereas the assay used full-length galectin-3. Reducing-end alpha and beta models are separate sensitivity conditions; solution ensemble weights are unknown and the six modeled inputs are not six independent assay ligands.

| Ligand | Heavy atoms | Measured KD ±SEM µM | Alpha mean Top1 energy kcal/mol | Beta mean Top1 energy kcal/mol |
| --- | --- | --- | --- | --- |
| LacNAc | 26 | 24.0 ± 3.6 | -5.873 | -6.011 |
| galactose | 12 | 3100.0 ± 1100.0 | -4.136 | -4.549 |
| lactose | 23 | 93.5 ± 9.1 | -5.036 | -5.251 |

Correlations compare higher-is-stronger scores (negative Vina energy, or heavy-atom count) with measured −log10(KD in M), separately by anomer. They describe three point estimates; no significance or generalization claim follows.

| Anomer | Score | Ligands | Spearman ρ | Kendall τ-b |
| --- | --- | --- | --- | --- |
| alpha | vina | 3 | 1.000 | 1.000 |
| alpha | heavy_atom_baseline | 3 | 1.000 | 1.000 |
| beta | vina | 3 | 1.000 | 1.000 |
| beta | heavy_atom_baseline | 3 | 1.000 | 1.000 |

Size alone ranks this panel perfectly, so the raw Vina rank does not demonstrate recognition beyond size. All compounds are neutral; a charge-only score is constant. Exact size/charge matching is unavailable in this three-ligand panel. No dissociation constant is inferred from docking energy.

## Execution resources

Wall and CPU values are sums of recorded final-attempt run times; wall sums are not elapsed study duration when jobs overlap. Prior-attempt counts are shown separately; their time is not included here. Each frozen run used two requested CPUs. Returned poses are retained candidates, not an internal sampling count.

| Case | Condition | Eligible runs | Completed runs | Returned poses | Wall sum s | CPU sum s | Prior attempts |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 154L | holo_complete | 3 | 3 | 60 | 542.631 | 942.082 | 0 |
| 1GNY | holo_complete | 3 | 3 | 60 | 588.931 | 1031.025 | 0 |
| 1KJL | holo_complete | 3 | 3 | 60 | 203.816 | 259.334 | 0 |
| 1PWB | holo_complete | 3 | 3 | 60 | 97.831 | 191.934 | 0 |
| 2J1T | holo_complete | 3 | 3 | 60 | 716.270 | 1271.881 | 0 |
| 2VUZ | holo_complete | 3 | 3 | 60 | 6385.620 | 12488.909 | 1 |
| AAV9 | apo_complete | 3 | 3 | 60 | 22.422 | 42.770 | 0 |
| AAV9 | apo_monomer | 3 | 3 | 59 | 23.065 | 43.485 | 0 |
| AAV9 | holo_complete | 3 | 3 | 60 | 22.918 | 43.044 | 0 |
| AAV9 | holo_monomer | 3 | 3 | 59 | 22.254 | 43.011 | 0 |
| gal3_LacNAc_alpha | holo_complete | 3 | 3 | 60 | 203.663 | 259.442 | 0 |
| gal3_LacNAc_beta | holo_complete | 3 | 3 | 60 | 149.725 | 247.725 | 0 |
| gal3_galactose_alpha | holo_complete | 3 | 3 | 60 | 43.587 | 50.313 | 0 |
| gal3_galactose_beta | holo_complete | 3 | 3 | 60 | 29.178 | 49.647 | 0 |
| gal3_lactose_alpha | holo_complete | 3 | 3 | 60 | 158.569 | 207.305 | 0 |
| gal3_lactose_beta | holo_complete | 3 | 3 | 60 | 162.606 | 210.729 | 0 |

Included totals: **48/48 completed runs**, **958 returned poses**, **9373.087 wall-seconds summed** and **17382.635 CPU-seconds summed**. Preparation failures remain in the seven-case table rather than disappearing from scope.

## Remaining scientific closure

The next evidence gates remain independently qualified structural references and assemblies, frozen family splits and training-overlap audits, blind/predicted-site methods and comparable information budgets, qualified tested alternatives, held-out execution and family-level uncertainty. The present results support implementation diagnostics and bounded failure observations only. They do not complete the original scientific question.

## Source provenance

SHA256 of the sorted path/hash manifest below: `97a2a2c41299432517ca135cfa800b1ff0c8a962b2fed81570c62b1220a6cae9`. Paths are repository-relative. All listed artifacts were read or verified during generation; no timestamp or machine-specific path is inserted into this report.

| Source | SHA256 |
| --- | --- |
| `aav_benchmark/report.py` | `b789461046d80994166778a2a94833c19ef5fd3bfc500186d6f9013b82e502c9` |
| `aav_benchmark/statistics.py` | `f8bb4e93865c14947c499493f9c8ebd4a55cb82c93d502317c3984e75202e64e` |
| `benchmarks/development.json` | `6c86cf008e9d1cbeea48fe847767abcdcae3b6410a475607e16f7c2bb63d4148` |
| `benchmarks/development_preparation.json` | `a843747a152d84627811974a418486aaa781047048aca73da1e1194b3ba1c6a2` |
| `benchmarks/external_run.json` | `0dc06423c0899339284d1e90c1ed07cd47787db5dd45c9099c512b91216bc355` |
| `benchmarks/gal3_run.json` | `0260057280030c207bc3c4f38795e1a142757a15702b6ee863714032f5f0dcd5` |
| `benchmarks/galectin3_itc_panel.json` | `ee3500ca9ab8fb838dc44a3c086068309155ee0fa9f0452a62c050b5d1219043` |
| `results/development/154l/configuration.json` | `6c86cf008e9d1cbeea48fe847767abcdcae3b6410a475607e16f7c2bb63d4148` |
| `results/development/154l/holo_complete_104729.pdbqt` | `eb3c98ba277077bc98c0ec8281d3fcebaf488bff2ac898316b3ed04e1b318020` |
| `results/development/154l/holo_complete_130363.pdbqt` | `662ad1a24cfe0b46ed75f0ad6e46c289bf65b85ce09d47039779179436711d26` |
| `results/development/154l/holo_complete_155921.pdbqt` | `bf92513d7975f18f840eb780f511c887b1fec3aed44211b294257d726a3aeb33` |
| `results/development/154l/pose_metrics.csv` | `b4fb545937dc55b66aa8b471ab6db392a76f059a809661d280b775107d479564` |
| `results/development/154l/runs.json` | `30aad06fd4a79b212162298db9083a42cce872a3672caca96a6a57c7c77c2ef0` |
| `results/development/154l/summary.json` | `01f6c09bd0c5277c8ca7e4560bf061e94ebd4f886771e37072ff00bc1bf237dc` |
| `results/development/1GNY/configuration.json` | `0dc06423c0899339284d1e90c1ed07cd47787db5dd45c9099c512b91216bc355` |
| `results/development/1GNY/holo_complete_104729.pdbqt` | `ab518bad0a6bc4eb0f6437722b3f25e8a453c67a2ecc104e7f2507c16a0490ab` |
| `results/development/1GNY/holo_complete_130363.pdbqt` | `29c024c4ef8d3e7bb283c15014821413b746670d8e33a9a782f5e0d77609f2c2` |
| `results/development/1GNY/holo_complete_155921.pdbqt` | `1ee591b836f5acf13c75b434428dfc478efa5b0d3467e38649b94f17cdb5da60` |
| `results/development/1GNY/pose_metrics.csv` | `8a9036e899a977923a7d254c1f2363f152ca1894d12bd406033b6dc4dee75792` |
| `results/development/1GNY/runs.json` | `13df7daa2751655493eca1d1c15b75a5cc21ef22952fe6d86404de303db4ca11` |
| `results/development/1GNY/summary.json` | `3ea2c213c3aa61f9eccb6086f4274518821b6a1b1e26267f9f1efaefb9ce631c` |
| `results/development/1KJL/configuration.json` | `0dc06423c0899339284d1e90c1ed07cd47787db5dd45c9099c512b91216bc355` |
| `results/development/1KJL/holo_complete_104729.pdbqt` | `58f32edf4cd5f1037ec3f9083ac6739d55941764a2c5ef8f64d5496102900f29` |
| `results/development/1KJL/holo_complete_130363.pdbqt` | `6f9bd7239caf541c0b9ae48283a08166293426eb76c41f7dffb5c3737d5728f5` |
| `results/development/1KJL/holo_complete_155921.pdbqt` | `ef3adf523f6f07926e2a5bf8c8b1121418d36a41851bc2e289c95ca12b5a34c3` |
| `results/development/1KJL/pose_metrics.csv` | `db7f4fbe021dce82b4a4805ed8100eb458a1b22bed0fa892e5e52eb5bf21efdb` |
| `results/development/1KJL/runs.json` | `df608bffa64134eae2431c459185b20bee4d28e6de1b8e5c87cc5ed313fd64a0` |
| `results/development/1KJL/summary.json` | `18c9079e9e821e4af9e6fd15033a0e78fe5f33fcbd8ce07ce40e9655001a08e2` |
| `results/development/1PWB/configuration.json` | `0dc06423c0899339284d1e90c1ed07cd47787db5dd45c9099c512b91216bc355` |
| `results/development/1PWB/holo_complete_104729.pdbqt` | `4ee5b7bec82d38cd322002d893ee457e84a53d5217853f370b6b7451431c7ce7` |
| `results/development/1PWB/holo_complete_130363.pdbqt` | `1d7a9c00bccca10ca2745a5ba7b771b6053e8095f8683456d87d23e5a0300986` |
| `results/development/1PWB/holo_complete_155921.pdbqt` | `91ddedb041ade6c709239a89127a1fb6df8fe1a3f07a71330381988d415f083b` |
| `results/development/1PWB/pose_metrics.csv` | `ebd388e1b31129e8386e923f71417b4158e5a9b6d5214f1b503c6b7b6148d84f` |
| `results/development/1PWB/runs.json` | `12162376e41988a2c795f9f4796bc08031451e932704e44939e2a263b389f605` |
| `results/development/1PWB/summary.json` | `8af70346a98ad51667b12e9f3959e67c8f01e5744c7e457cd4589f6eb929233a` |
| `results/development/2J1T/configuration.json` | `0dc06423c0899339284d1e90c1ed07cd47787db5dd45c9099c512b91216bc355` |
| `results/development/2J1T/holo_complete_104729.pdbqt` | `d670493c5e21d8c2b19aebc5dbae18d847b7feb6277a4788aecb7482200aa5a2` |
| `results/development/2J1T/holo_complete_130363.pdbqt` | `2b9c9ff064fa1ee66d6a25bcd5d9cfddd91939c952a88b54a1a1e58ba46877ab` |
| `results/development/2J1T/holo_complete_155921.pdbqt` | `4d29acb6c9b9f7f54fd04872072a9842b6a0e2a9691e65e5761200510a700ebd` |
| `results/development/2J1T/pose_metrics.csv` | `afe2f02e3232b64d6e772a79793f29ea6868b9d8d6ab572ad73d51701f72c103` |
| `results/development/2J1T/runs.json` | `35de4f6b0552f767a9b4bdc5bfe08a30c5a1e7c38b246c1e6e6b300d9ff723a0` |
| `results/development/2J1T/summary.json` | `9fc8ab3fc8b7bc38db43efd5674286025fec9f5f362d5ebaa5621422ad7a495e` |
| `results/development/2VUZ/configuration.json` | `0dc06423c0899339284d1e90c1ed07cd47787db5dd45c9099c512b91216bc355` |
| `results/development/2VUZ/holo_complete_104729.pdbqt` | `61ffd897b87f58a82aab9a1e0fe534d3b866eecc865805c6236d788e3698e921` |
| `results/development/2VUZ/holo_complete_130363.pdbqt` | `3be6808f81b53b6101ee395d80b1e654d3c770605e07d4a62cbe96d3f4fa5dbc` |
| `results/development/2VUZ/holo_complete_155921.pdbqt` | `1c92e099a81c56ed95d9848deb68e49efae79ba56ae65c879214bce474607179` |
| `results/development/2VUZ/pose_metrics.csv` | `cab8f042d21cb863fbd7ef49d1508002af91db56366fc53d744eec6721485a1a` |
| `results/development/2VUZ/runs.json` | `169083d3baa181f8c77c18ca6f61b7ea27da2239aaba93de02fcab555ebb13e2` |
| `results/development/2VUZ/summary.json` | `d43ed3f3a0daaa538a79600a93e4ffe5720a01f87c98afbd045da5e634528200` |
| `results/development/aav9/apo_complete_104729.pdbqt` | `cfafc70ba94c7b90880a96615358523d5a84d6b92d8abbc0217d9b024f4ac5a8` |
| `results/development/aav9/apo_complete_130363.pdbqt` | `1b40a2e2aecc4fce16f3b79c0a824b503d3473d425626bcffd863791baaf12ff` |
| `results/development/aav9/apo_complete_155921.pdbqt` | `eef852cc19fece5b3caef7e81f1def140dd81860dc4c5ee137075cf5b3276a25` |
| `results/development/aav9/apo_monomer_104729.pdbqt` | `8a3b2fa2e77dd7569d1d5a1a45a1c2e71424ac527dcc0883db6ac5c8b362f560` |
| `results/development/aav9/apo_monomer_130363.pdbqt` | `f4836af143cfe121566256e4abaeef375ad6b759a83b55984f37a6593ca131fe` |
| `results/development/aav9/apo_monomer_155921.pdbqt` | `1ddaf2140d91032863b2860c1723a6e0ac4aa79ccb73d8dfa6b00608631851bf` |
| `results/development/aav9/configuration.json` | `6c86cf008e9d1cbeea48fe847767abcdcae3b6410a475607e16f7c2bb63d4148` |
| `results/development/aav9/holo_complete_104729.pdbqt` | `2e2c5d7857f757f5ead40f44bd7b859a7711ae0bca8021c54e7dfbef6d956031` |
| `results/development/aav9/holo_complete_130363.pdbqt` | `87a2645bc2e2bbfa01e598743b642bb5fbbb75599fc322f57f843fc89c1021fa` |
| `results/development/aav9/holo_complete_155921.pdbqt` | `544e9b2e45c1529866edba34693446cd5577eb3eea0d4f5670bc95ae670e0347` |
| `results/development/aav9/holo_monomer_104729.pdbqt` | `dcda3952dda15a8a18c0af978fdc25e2a0911a29df32700608b365a5ed96f3b0` |
| `results/development/aav9/holo_monomer_130363.pdbqt` | `c83a2274b4bf59dd5583b2133c6fa2c961454d62168cacd11523256f85125a85` |
| `results/development/aav9/holo_monomer_155921.pdbqt` | `ff878d026ffb6bfcf4b8196c111f9ee2540fb67be64b45c47aef216f713f1f7a` |
| `results/development/aav9/native_control/apo_complete_score_only.log` | `0a17fd0d8bf5fb0437c764d965d8f5f84fe029cf148d481fc66e96da67905a21` |
| `results/development/aav9/native_control/apo_monomer_score_only.log` | `9e6aa0a27695fe781ef8f9a5434255b8f0452e2d2262961357cca2fbf4a54f2b` |
| `results/development/aav9/native_control/holo_complete_score_only.log` | `868865e7928ddc042c1ae9e59916f96d1d06717d6a7f3924467f797b326b1ac4` |
| `results/development/aav9/native_control/holo_monomer_score_only.log` | `3abb9d60eb5124c67231a0cb52185a1479991d86ea23c04aedef090afa848cdd` |
| `results/development/aav9/native_control/ligand_preparation.json` | `ce8edfaa7d9ab215ccd6a0263f58c1cd79035e563aa1236882fb973fb6c126fe` |
| `results/development/aav9/native_control/native_ligand.pdbqt` | `e4c4678a8c9d67f388c9c5a5d9a7e6c84ceac9c6fd051c8cc2e107bb73cbd5af` |
| `results/development/aav9/native_control/native_ligand.sdf` | `e1d75b9b309c5e144ef5fa4494067d4d7605f484649914b5cccc48def4c8b58f` |
| `results/development/aav9/native_control/summary.json` | `8fd273a00d198afc95cd20c9d7423925c258e1c0a0e38583bc5dd13232287274` |
| `results/development/aav9/pose_metrics.csv` | `601dea028e5621e220ad41c84794061097dc8c54341dd15f4033f4e1cd402766` |
| `results/development/aav9/runs.json` | `f2368cc9f0afd8d82e2df1ca6998fff39df1b18918134e6c1b80e1ef45fd607e` |
| `results/development/aav9/summary.json` | `a41b27583bd4b17c8187dc3d11c701c5668ceca335291a98ffb2f870450bdf35` |
| `results/development/gal3_LacNAc_alpha/configuration.json` | `0260057280030c207bc3c4f38795e1a142757a15702b6ee863714032f5f0dcd5` |
| `results/development/gal3_LacNAc_alpha/holo_complete_104729.pdbqt` | `0a5482ce87c635496ebd6222ccb1c7ccc6fdd794a02731b17be82248b9eee6c3` |
| `results/development/gal3_LacNAc_alpha/holo_complete_130363.pdbqt` | `56945157375c746f89beb22377d2ec0c14bdea7e61db4d928c1bd32a41e7c3c2` |
| `results/development/gal3_LacNAc_alpha/holo_complete_155921.pdbqt` | `9b0d30f8c729ac1e0d91460982869eca0e27a96a29587a94521fe605206b3a76` |
| `results/development/gal3_LacNAc_alpha/runs.json` | `20ec951317696ac78f75f97bbc22fc0e81752a4d09836b01d4b52487fed0c99b` |
| `results/development/gal3_LacNAc_beta/configuration.json` | `0260057280030c207bc3c4f38795e1a142757a15702b6ee863714032f5f0dcd5` |
| `results/development/gal3_LacNAc_beta/holo_complete_104729.pdbqt` | `58f32edf4cd5f1037ec3f9083ac6739d55941764a2c5ef8f64d5496102900f29` |
| `results/development/gal3_LacNAc_beta/holo_complete_130363.pdbqt` | `6f9bd7239caf541c0b9ae48283a08166293426eb76c41f7dffb5c3737d5728f5` |
| `results/development/gal3_LacNAc_beta/holo_complete_155921.pdbqt` | `ef3adf523f6f07926e2a5bf8c8b1121418d36a41851bc2e289c95ca12b5a34c3` |
| `results/development/gal3_LacNAc_beta/runs.json` | `164e5084f1cfdf827708f43008cc86334cb30d182dd5e88e9ec03f25c2590e3a` |
| `results/development/gal3_galactose_alpha/configuration.json` | `0260057280030c207bc3c4f38795e1a142757a15702b6ee863714032f5f0dcd5` |
| `results/development/gal3_galactose_alpha/holo_complete_104729.pdbqt` | `c5871dc702014d75eb6927092b684d2dcdd78873a00de3a21282f696fa8c47d4` |
| `results/development/gal3_galactose_alpha/holo_complete_130363.pdbqt` | `d78955a2e5e0a3982d198085881ce7be1bd66e262c5d315b8a17cc29b1762c33` |
| `results/development/gal3_galactose_alpha/holo_complete_155921.pdbqt` | `e5eb99491533da28134104abce2c33edd4b714a31c678e2e63450e0cb734ded0` |
| `results/development/gal3_galactose_alpha/runs.json` | `fa1441b80bdaddb50804933b708cd0b669771137cd9e992b4fa6d4904350320b` |
| `results/development/gal3_galactose_beta/configuration.json` | `0260057280030c207bc3c4f38795e1a142757a15702b6ee863714032f5f0dcd5` |
| `results/development/gal3_galactose_beta/holo_complete_104729.pdbqt` | `a036b135650c47b159b0a69787364b30a80c0979807a2da54f66979c8db30382` |
| `results/development/gal3_galactose_beta/holo_complete_130363.pdbqt` | `d5324fccbe96a8b3cb47b8569a672ea61e035345584e631d7ec562340e18868f` |
| `results/development/gal3_galactose_beta/holo_complete_155921.pdbqt` | `a78a741a62ce40c64833d97ae279dc12935e6cedf6f3aadbc3d4a57ce869e01b` |
| `results/development/gal3_galactose_beta/runs.json` | `d1282655fcba3d1d0e519444e7cfa102a10ef2ea49570129dc98cb3447d65b79` |
| `results/development/gal3_lactose_alpha/configuration.json` | `0260057280030c207bc3c4f38795e1a142757a15702b6ee863714032f5f0dcd5` |
| `results/development/gal3_lactose_alpha/holo_complete_104729.pdbqt` | `c154bccc2b874584edb8e64c88f35e5303be2e6c39267c67eff2a0ce688e6508` |
| `results/development/gal3_lactose_alpha/holo_complete_130363.pdbqt` | `d09e083eb4d80e65e3315c1eb3802ac01b7b24fc327144cd5873fdbe3f6ef196` |
| `results/development/gal3_lactose_alpha/holo_complete_155921.pdbqt` | `c02f4b5d078b82373c4912909c1aa1664e80ec2b04cf49bcf1f0850f1d8f0e5c` |
| `results/development/gal3_lactose_alpha/runs.json` | `d2e6f46283446f0d2c2840f0b2a8baf3bf09ad4c5e9f4d67d3f012167a7b7279` |
| `results/development/gal3_lactose_beta/configuration.json` | `0260057280030c207bc3c4f38795e1a142757a15702b6ee863714032f5f0dcd5` |
| `results/development/gal3_lactose_beta/holo_complete_104729.pdbqt` | `60eb3cf204f1486a8140e0999773ec1605ffd55c3671869c4f12ac94316acfb2` |
| `results/development/gal3_lactose_beta/holo_complete_130363.pdbqt` | `ed1a5c5be6638ec3c9c38c5f69bf208438a13d14c9ecbae145eb105cf222cf11` |
| `results/development/gal3_lactose_beta/holo_complete_155921.pdbqt` | `ab69511433998f7b94b8e6bfcffd9277ff64e79bb8c4d15ed4e65125865f6bf7` |
| `results/development/gal3_lactose_beta/runs.json` | `0268d34e036f8f87c86491d18b81c055c527e815e21f631dc605f69a8dfc24ef` |
| `results/development/gal3_panel_summary.json` | `908d962fdd081f581caded4cb894d0ce61e3d56d41c78b7b430f42df01615783` |
| `scripts/regenerate_report.py` | `406f199a286a33f5cb6801d17e7129e6ea720f9338268c9830fecfeca8feacde` |
