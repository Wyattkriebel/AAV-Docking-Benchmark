# POL: an independent natural epimer contrast

This frozen known-site analysis compares measured D-galactose binding with actually tested D-glucose showing no detectable ITC binding. NB is categorical; no numerical Kd bound is assigned. This is one contrast in the POL/CBM67 family.

The two sugars differ only at C4 and share formula C6H12O6, 12 heavy atoms and neutral model charge.

Raw verification: **96/96 successful runs**; 0 failed, 0 running and 0 missing.

**Primary N-site, holo, calcium-retained condition: D-galactose scores lower.** Equal-anomer means are -4.2045 kcal/mol for D-galactose and -4.1772 kcal/mol for D-glucose. The glucose-minus-galactose difference is 0.0273 kcal/mol.

Seeds are averaged within each modeled anomer, then the two anomers receive equal weight. This is a computational heuristic, not measured population weighting. Furanose, open-chain and solution anomer fractions are unknown. Sites, receptor states, calcium profiles, anomers and seeds are sensitivity conditions rather than independent biological replicates.

The ranking direction reverses in **6 of 7 prespecified sensitivity conditions**. These conditions are descriptive checks of this one contrast, not independent biological replications.

The matched alpha and beta model pairs also favor different sugars in the primary condition. The equal-weight primary ranking is sensitive to this unmeasured solution-mixture choice.

## Separate conditions

| State | Site | Calcium profile | Gal mean | Glc mean | Glc-Gal | Gal min/max anomer | Glc min/max anomer |
| --- | --- | --- | --- | --- | --- | --- | --- |
| holo | N | Ca_retained | -4.2045 | -4.1772 | 0.0273 | -4.2200 / -4.1890 | -4.2063 / -4.1480 |
| holo | N | dry | -3.9450 | -4.0857 | -0.1407 | -3.9617 / -3.9283 | -4.1640 / -4.0073 |
| holo | C | Ca_retained | -4.4187 | -4.4405 | -0.0218 | -4.4253 / -4.4120 | -4.4837 / -4.3973 |
| holo | C | dry | -4.1518 | -4.1733 | -0.0215 | -4.2030 / -4.1007 | -4.2220 / -4.1247 |
| apo | N | Ca_retained | -3.7525 | -3.8322 | -0.0797 | -3.7700 / -3.7350 | -3.9243 / -3.7400 |
| apo | N | dry | -3.5635 | -3.6362 | -0.0727 | -3.6257 / -3.5013 | -3.7577 / -3.5147 |
| apo | C | Ca_retained | -4.0940 | -4.2843 | -0.1903 | -4.1270 / -4.0610 | -4.3387 / -4.2300 |
| apo | C | dry | -3.9645 | -3.8838 | 0.0807 | -4.0380 / -3.8910 | -3.8970 / -3.8707 |

All energies are kcal/mol; lower scores are favorable. `dry` removes calcium; both profiles omit explicit waters.

## Every anomer and seed

| State | Site | Profile | Species | Anomer | Seeds 104729 / 130363 / 155921 | Mean |
| --- | --- | --- | --- | --- | --- | --- |
| holo | N | Ca_retained | D_galactose | alpha | -4.190, -4.189, -4.188 | -4.1890 |
| holo | N | Ca_retained | D_galactose | beta | -4.223, -4.229, -4.208 | -4.2200 |
| holo | N | Ca_retained | D_glucose | alpha | -4.133, -4.236, -4.250 | -4.2063 |
| holo | N | Ca_retained | D_glucose | beta | -4.122, -4.162, -4.160 | -4.1480 |
| holo | N | dry | D_galactose | alpha | -3.928, -3.927, -3.930 | -3.9283 |
| holo | N | dry | D_galactose | beta | -3.962, -3.961, -3.962 | -3.9617 |
| holo | N | dry | D_glucose | alpha | -4.014, -4.006, -4.002 | -4.0073 |
| holo | N | dry | D_glucose | beta | -4.159, -4.169, -4.164 | -4.1640 |
| holo | C | Ca_retained | D_galactose | alpha | -4.433, -4.383, -4.420 | -4.4120 |
| holo | C | Ca_retained | D_galactose | beta | -4.424, -4.425, -4.427 | -4.4253 |
| holo | C | Ca_retained | D_glucose | alpha | -4.393, -4.390, -4.409 | -4.3973 |
| holo | C | Ca_retained | D_glucose | beta | -4.457, -4.544, -4.450 | -4.4837 |
| holo | C | dry | D_galactose | alpha | -4.178, -4.156, -4.275 | -4.2030 |
| holo | C | dry | D_galactose | beta | -4.104, -4.097, -4.101 | -4.1007 |
| holo | C | dry | D_glucose | alpha | -4.129, -4.120, -4.125 | -4.1247 |
| holo | C | dry | D_glucose | beta | -4.220, -4.220, -4.226 | -4.2220 |
| apo | N | Ca_retained | D_galactose | alpha | -3.732, -3.732, -3.741 | -3.7350 |
| apo | N | Ca_retained | D_galactose | beta | -3.763, -3.777, -3.770 | -3.7700 |
| apo | N | Ca_retained | D_glucose | alpha | -3.740, -3.733, -3.747 | -3.7400 |
| apo | N | Ca_retained | D_glucose | beta | -3.891, -3.950, -3.932 | -3.9243 |
| apo | N | dry | D_galactose | alpha | -3.513, -3.476, -3.515 | -3.5013 |
| apo | N | dry | D_galactose | beta | -3.633, -3.621, -3.623 | -3.6257 |
| apo | N | dry | D_glucose | alpha | -3.532, -3.497, -3.515 | -3.5147 |
| apo | N | dry | D_glucose | beta | -3.760, -3.770, -3.743 | -3.7577 |
| apo | C | Ca_retained | D_galactose | alpha | -4.064, -4.058, -4.061 | -4.0610 |
| apo | C | Ca_retained | D_galactose | beta | -4.124, -4.128, -4.129 | -4.1270 |
| apo | C | Ca_retained | D_glucose | alpha | -4.220, -4.237, -4.233 | -4.2300 |
| apo | C | Ca_retained | D_glucose | beta | -4.351, -4.333, -4.332 | -4.3387 |
| apo | C | dry | D_galactose | alpha | -3.888, -3.892, -3.893 | -3.8910 |
| apo | C | dry | D_galactose | beta | -4.051, -4.017, -4.046 | -4.0380 |
| apo | C | dry | D_glucose | alpha | -3.872, -3.867, -3.873 | -3.8707 |
| apo | C | dry | D_glucose | beta | -3.904, -3.894, -3.893 | -3.8970 |

## Verification and scope

The audit checks frozen source/configuration/code snapshots, every preparation archive, exact CCD/SDF/PDBQT stereochemistry and atom maps, native dimer/calcium identities, each command, all raw pose hashes, every returned pose's heavy elements and stereochemistry, and raw log score tables. Recorded totals are 1919 returned poses, 1021.4 summed run wall-seconds and 1483.8 CPU-seconds.

The logs contain 1920 modes; 1 are outside the saved energy window. Vina prints all search-result modes before its PDBQT writer applies the frozen 10 kcal/mol window. The audit requires the saved poses to be the corresponding rank prefix, allowing only documented print rounding. [Vina 1.2.7 source](https://github.com/ccsb-scripps/AutoDock-Vina/blob/v1.2.7/src/lib/vina.cpp#L557-L613).

The holo A/B model preserves 5,034 observed protein heavy atoms; the independently solved apo dimer preserves 5,058 after coherent highest-occupancy Ser218B selection. Calcium-retained profiles add four native ions. A 2,704-atom protein-backbone alignment transfers centers into the apo frame. Native apo glycosylation, missing N-terminal residues and assay/construct correspondence limitations remain disclosed.

No AUROC/AUPRC or calibrated affinity is reported from one epimer pair. Fixed +2 Vina calcium approximates the native ions and does not establish accurate coordination energetics. Input coordinates are independently generated and supply no experimental pose truth.

Assay source: [Perduca et al. 2020, Figure 4 and ITC methods](https://doi.org/10.1093/glycob/cwaa006). [Source qualification](INDEPENDENT_SPECIFICITY_CANDIDATES.md), [frozen experiment](../benchmarks/pol_experiment.json), and [machine-readable evidence](../results/pol_specificity/summary.json).
