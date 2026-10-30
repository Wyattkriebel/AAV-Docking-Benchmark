# Additional independent evidence

Research checkpoint: 29 October 2026. A six-case Gal3 comparison can now run as a
literal reproduction of published classification labels. The recovered evidence
does not yet qualify an independently reconstructed biological specificity panel.
The frozen three-binder ITC panel is unchanged. Exact cells, chemical graphs,
checksums and unresolved qualifications are recorded in
`benchmarks/additional_assay_candidates.json`.

## Recovered Gal3 measurements and modeled graphs

[Carpenter et al. 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12680057/)
provide experimental LiGA measurements and a published model-comparison workbook.
Data S1 and the supplementary information were recovered through the
[official EuropePMC archive](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12680057/supplementaryFiles);
all archive members passed their CRC checks. The source workbook contains:

| Glycoside | Published “LiGA” field | Published label | Row |
|---|---:|---:|---:|
| α-L-Glc | 2.2520295213 | 1 | 12 |
| β-L-Glc | 1.0804303437 | 0 | 23 |
| β-D-Glc | 1 | 0 | 34 |
| α-L-Man | 0.3906077032 | 0 | 45 |
| α-D-Man | 0.913292175 | 0 | 56 |
| α-L-Rha | −0.5686447951 | 0 | 67 |

Values occupy columns E/F of `MCNet vs 4 other models.xlsx`, Sheet1. A second
workbook supplies source SMILES. All six parse with stereochemistry. The first
five modeled compounds each have 22 heavy atoms and net charge zero; α-L-Rha has
21 heavy atoms. This permits an exactly size/charge-matched source-reproduction
comparison. Size alone gives AUROC 0.6 across six source labels and 0.5 among the
five matched compounds. These are reference baselines, not new docking results.

Several source conflicts prevent stronger biological inference:

- The article broadly describes no Gal3 binding to the L-monomers, while the
  comparison workbook labels α-L-Glc positive.
- The supplementary definition describes a normalized output/input ratio, but
  one published LiGA value is negative. Its transformation or aggregation has
  not been independently reconstructed.
- Gal3 lectin-array measurements are absent; no second assay adjudicates these
  entries. The source's general OR rule must not be mistaken for two measured
  platforms in this particular subset.
- The source SMILES represent modeled spacer-bearing glycans. Actual experiments
  use multivalent glycophages, with conjugation chemistry and density that the
  modeled graphs do not encode.

The original Gal3 output-count datasets
[MOG](https://48hd.cloud/file/20157) and
[MOY](https://48hd.cloud/file/21789) were also recovered and hashed. They contain
replicate counts, but matching naive inputs, chemical barcode/density maps and
aggregation provenance remain missing from the recovered assets. Their stated
noncommercial distribution terms are recorded; the files are not copied into
this MIT-licensed repository. A new run can report fidelity to the six published
labels while keeping these limitations explicit. It cannot establish general or
free-sugar binding specificity.

## Direct negative control with unresolved chemical identity

[Stowell et al. 2008](https://pmc.ncbi.nlm.nih.gov/articles/PMC2442294/)
actually tested a biotinylated arabinose reference surface alongside lactose and
poly-LacNAc by SPR. The source calls arabinose nonbinding and reports 0.1–100 µM
galectin, 25°C and PBS pH 7 with 0.005% Tween20. Its control's D/L configuration,
ring/anomer and complete linker remain unidentified. Labeling an arbitrary free
arabinose would invent a chemical correspondence.

The same study reports Gal3 binding lactose by SPR while failing to recognize
immobilized lactose/LacNAc on its glycan array. Thus a low array signal cannot be
transferred to a free-ligand assay without preserving presentation.

## Conditional inhibition is a separate endpoint

[The SP-D/A2M study](https://pmc.ncbi.nlm.nih.gov/articles/PMC2859506/)
actually tested D-fucose at 25 mM as a poor-binding control alongside mannose and
maltose. The result supports conditional non-inhibition, not zero binding or a
numerical affinity bound. Its Methods and Figure 5A also disagree on protein and
coating concentrations; both statements are retained in the candidate ledger.

[McCartney et al. 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1450244/)
report CBM15 association constants of 20000 and 2000 M⁻¹ for xylohexaose and
cellohexaose, attributing the data to the original structural paper. Both bind;
the weaker alternative must not be assigned a negative label. Original assay
details are still needed before a new quantitative comparison.

## Local structural quality remains an explicit gate

The archived diagnostics for 2J1T show ligand RSCC 0.899–0.973, occupancy 1 and
RSR 0.049–0.091, but also 11 angle and two torsion outlier flags. 1GNY includes
terminal RSCC 0.808, occupancy 0.95 and geometric flags. Neither is free of caveats
on RSCC alone. An atomic-reference decision must inspect the flagged atoms,
complete graph/assembly and local 2Fo−Fc/difference density, with any accepted core
defined before evaluation. Current reference eligibility remains unchanged.
