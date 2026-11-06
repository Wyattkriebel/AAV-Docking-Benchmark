# A tested, size-matched specificity endpoint

Research checkpoint: 30 October 2026. A same-assay conditional specificity
panel is now available from
[Mahanti et al., DOI 10.1021/acsmedchemlett.9b00396](https://pmc.ncbi.nlm.nih.gov/articles/PMC6956361/).
The conservative cohort contains **six measured Gal9C ligands and six actually
tested C3 epimers reported as NB at 2 mM**. Every positive/negative pair has the
same formula, heavy-atom count, formal charge, and chemical connectivity; only
the sugar C3 configuration differs. This supplies the tested alternative class
that the three-binder Gal3 ITC panel lacks.

The endpoint is recognition in the source competitive fluorescence assay. NB
means the source found no binding at its highest tested concentration, 2 mM.
It does not mean zero affinity at every concentration, and it does not establish
a numerical `Kd > 2 mM` bound. No docking result or calibrated affinity claim is
made here. The frozen ITC panel and candidate ledger are unchanged.

The complete record, exact stereochemical SMILES, checksum provenance, source
conflicts, and a proposed computational analysis are in
[`benchmarks/specificity_resolution.json`](../benchmarks/specificity_resolution.json).

## Primary conservative cohort

Table 1 reports human galectin-9 C-terminal domain (Gal9C) affinities as means
and SD of 4–8 single–double point measurements. These are soluble methyl
β-D-pyranosides with fixed anomers, rather than free reducing sugars or phage
display constructs. The galacto compounds were previously synthesized, but the
paper explicitly evaluates both series together in this assay. Reference sugars
with affinity values imported from other publications are excluded.

| N-aryl substituent | Galacto compound | Gal9C Kd ± SD, µM | Tested gulo epimer | Gal9C result | Heavy atoms, each |
|---|---:|---:|---:|---|---:|
| phenyl | 1 | 678 ± 30 | 14 | NB at 2 mM | 19 |
| 4-methylphenyl | 4 | 770 ± 25 | 17 | NB at 2 mM | 20 |
| 3-methylphenyl | 5 | 670 ± 23 | 18 | NB at 2 mM | 20 |
| 3,4-dimethylphenyl | 6 | 130 ± 9 | 19 | NB at 2 mM | 21 |
| 3-chlorophenyl | 7 | 160 ± 6 | 20 | NB at 2 mM | 20 |
| 3,4-difluorophenyl | 8 | 220 ± 7 | 21 | NB at 2 mM | 21 |

Use the source's finite-affinity versus NB distinction. An untested compound or
an ND entry is never a negative. There are no ND entries in the Gal9C cohort.
The source's [Chart 1 and Scheme 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC6956361/)
and the [official supporting information](https://acs.figshare.com/doi/suppl/10.1021/acsmedchemlett.9b00396.s001)
identify the graphs. The SI reports >95% analytical HPLC purity for tested
compounds.

RDKit checks all 16 transcribed graphs, including the optional methoxy pairs:
all expected formulas and heavy-atom counts match, every molecule has five
defined stereocenters, and every paired graph differs at exactly one center,
C3. Removing stereochemistry makes each pair identical. Both the conservative
12-row cohort and the expanded 16-row cohort give a heavy-atom-count AUROC of
0.5 in either orientation; all graph formal charges are zero.

## Assay and construct provenance

The 2019 SI S2 explicitly adopts earlier galectin-specific conditions. Its
reference 3,
[Delaine et al. 2016](https://www.pure.ed.ac.uk/ws/portalfiles/portal/27411623/Delaine_et_al_2016_ChemBioChem.pdf),
specifies PBS, 20°C, Gal9C at 2.0 µM, and a fluorescein thiodigalactoside probe at
0.10 µM, with 485/520 nm excitation/emission. The probe's complete name is
retained in the JSON. The published version corresponding to SI reference 4 is
[Pal et al. 2019, DOI 10.3390/ijms20153786](https://pmc.ncbi.nlm.nih.gov/articles/PMC6696278/);
its section 3.20 also retains these earlier galectin-specific conditions.

The 2016 method encodes human Gal9C residues 205–355 with thioredoxin in pET-32
Ek/LIC. The recovered method does not describe cleavage of that fusion. The
2019 paper does not repeat an assayed sequence/tag or batch record. Thus the
residue range is method provenance, while exact 2019 construct identity remains
an approximation to disclose. PBS pH and ionic composition are not specified.
These gaps do not turn NB into an untested label, but they limit experimental
reconstruction.

The source uses 3NV3 for Gal9C molecular dynamics. It is a candidate structure
for a later receptor decision; that use does not establish atomic-reference
quality or a site-blind information regime.

## Source conflicts are explicit

The conservative cohort omits compounds 2/15 and 3/16. Main Scheme 1 identifies
15 as the meta-methoxy guloside and 16 as the para-methoxy guloside, while the
SI synthesis heading reverses those names and contains additional numbering
errors. Both separated gulo regioisomers were tested, and both are NB for
Gal9C. Consequently the two-graph Gal9C class set is invariant to the swap, but
the expanded 16-compound cohort is retained only as a named sensitivity check.
Gal9N numerical measurements must not be transferred through those ambiguous
aliases.

Main molecular-dynamics prose calls 6 and 18 corresponding epimers. Main Chart
1 makes 6 dimethyl, while Scheme 1 makes 18 monomethyl and 19 dimethyl. The SI
compound-19 heading, two methyl NMR signals, and formula `C15H23NO5` agree with
Scheme 1; its synthesis body also mistakenly says “afford 20.” The chemically
matched pair is **6/19**. No claim about source MD results is transferred from
18 to 19.

## Analysis that can now be run

Freeze the receptor, construct approximation, site-information regime, sugar
ring generation, protonation, random seeds, and score aggregation before new
docking. Preserve the five stereocenters and fixed methyl β-anomer. Use all six
conservative epimer pairs. Report individual scores and how often the measured
galacto ligand scores more favorably than its tested gulo counterpart, with a
tie rule declared in advance. Report descriptive AUROC alongside the size and
charge baselines; any score threshold must be chosen independently of this
panel. Keep NB categorical and never substitute a numerical Kd.

This is one domain and one synthetic glycoside series selected after inspecting
the literature. It supports an exploratory, chemically controlled specificity
test. It does not supply held-out receptor-family generalization or a
confirmatory benchmark by itself.

## Gal3 fallback and Carpenter audit

[Zetterberg et al. 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9574852/), Table 1,
provides a second actual tested contrast: compounds 3 and 4 have finite FP Kd
values, while 5a and 5b are each right-censored at `>100 µM`. The official CSV
swaps the 5a/5b structures. Main synthesis names, NMR, HRMS, and SI headings
identify 5a as penta-O-methyl and 5b as hexa-O-methyl; both have the same censoring,
so the correction changes no affinity class. Exact corrected graphs and evidence
are retained in the JSON. This fallback is size-confounded: its positives have
44/49 heavy atoms and negatives 50/51. Smaller size alone gives AUROC 1.0. The
paper also omits a complete explicit FP recipe and construct record.

For [Carpenter et al. 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12680057/),
both named naïve-input files were recovered and CRC-checked. This closes the
specific missing-input retrieval gap in the earlier research checkpoint.

| Role | Public file | Count columns | Data rows excluding total |
|---|---|---:|---:|
| MOG input | [20250131-87MOGooNA-IF](https://48hd.cloud/file/20163) | 11 | 1,295,021 |
| MOG Gal3 output | [20250131-87MOGgaPANN-IF](https://48hd.cloud/file/20157) | 3 | 151,125 |
| MOY input | [20250716-87MOYooNA-EG](https://48hd.cloud/file/21644) | 4 | 206,198 |
| MOY Gal3 output | [20250718-87MOYgaPANN-IF](https://48hd.cloud/file/21789) | 3 | 176,935 |

The tables identify `MOG0.SDB0`–`SDB35` and `MOY0.SDB0`–`SDB35`, with many
unassigned sequence rows, but do not supply chemical graph/density definitions.
The choice of input columns and the aggregation that produced the six comparison
workbook values remain undocumented in the recovered material. Complete
recursive inventories of the authors' MCNet and mirror-image-glycans repositories
were checked. None of the six MCNet processing scripts contains the searched
LiGA, MOG/MOY, barcode, edgeR, or logFC terms. Checksums and repository revisions
are recorded; this is an inventory result, not proof that data exist nowhere
else.

The negative α-L-Rha value still does not establish a transformation, and the
α-L-Glc positive workbook label still conflicts with the body/Figure S15
generalization. No biological Carpenter labels were promoted. The newly
recovered 48HD tables remain outside this MIT repository under their source
terms; their hashes, dimensions, and public retrieval pages provide the audit.
