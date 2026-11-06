# Independent specificity sources

Source inventory and POL qualification frozen **3 November 2026**, before any
POL docking or score inspection. This adds one independent protein family with
an actually tested natural epimer alternative. It does not establish general
reliability across families. Existing assay records are unchanged.

The complete qualification record is
[`independent_specificity_candidates.json`](../benchmarks/independent_specificity_candidates.json).
The computational experiment must have its own freeze.

## Qualified source contrast: oyster mushroom lectin POL

[Perduca et al. 2020](https://doi.org/10.1093/glycob/cwaa006), Figure 4A-C and
p.557, measures free galactose and actually titrates its C4 epimer glucose.
The D configurations are shown in Figure 4E-F. These are soluble reducing
sugars, with no carrier or linker.

| Source species | ITC observation | Primary role |
|---|---|---|
| D-galactose | Binding; main text gives approximately 0.78 mM Kd | Measured binder |
| D-glucose | No binding detected in the reported titration | Tested alternative |
| Melibiose | Binding; approximately the same Kd as galactose | Retained source inventory; larger-ligand sensitivity only |

The method on p.560 uses 190 µL of 100 µM recombinant POL in 20 mM Tris/HCl,
pH 7.5, with 0.1 mM CaCl2 at 25°C. Twenty 2.5 µL injections contain 10 mM
sugar in the same buffer; stirring is 250 rpm. Ligand-into-buffer dilution
heats are subtracted. **10 mM is the titrant concentration**, not a claimed
terminal cell concentration or Kd bound. No numerical detection threshold is
given. NaCl is absent from this ITC recipe, while prior gel filtration used
0.15 M NaCl; the residual salt/exchange history is unspecified.

Separate STD-NMR experiments corroborate galactose recognition and lack of a
glucose signal. Their conditions differ from ITC and remain separate evidence.
The official page lists `SI_cwaa006`; its bytes and Table SIII were not
recovered. Exact uncertainties and replicate counts are therefore unassigned.

## Exact model graphs and reducing-sugar limitation

The archived [CCD definitions](../data/archive/independent_specificity_sources/manifest.json)
provide these explicit model states. Each has formula `C6H12O6`, 12 heavy atoms,
graph formal charge zero, and five defined stereocenters, verified with RDKit.
The matched-anomer galactose/glucose pairs differ only at C4.

| State | CCD | Isomeric SMILES |
|---|---|---|
| α-D-galactopyranose | GLA | `OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O` |
| β-D-galactopyranose | GAL | `OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O` |
| α-D-glucopyranose | GLC | `OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O` |
| β-D-glucopyranose | BGC | `OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O` |

Four model states represent **two experimental species**. The source does not
measure their anomer populations or establish locked anomers. Furanose and
open-chain fractions are unknown. An equal-anomer score mean is a declared
model heuristic; it is not population weighting. Best/worst-anomer sensitivity
must accompany it. No independent experimental label is assigned to each model
state.

## Receptor correspondence and the preselected known site

[6T1D](https://www.rcsb.org/structure/6T1D) is the same study's recombinant
melibiose complex. Assembly 1 is the A/B disulfide-linked dimer, with two sites
per chain and four per dimer. The asymmetric unit contains three dimers. Sites
and chains do not become independent families.

The deposited sequence contains 356 residues: the annotated `GSHM` expression
remnant followed by UniProt E7E2M2 positions 22-373. Its 352-residue core is
identical to the corresponding core in native apo
[6T0Q](https://www.rcsb.org/structure/6T0Q). Both entries annotate D31G, S43R,
N51S, and A52T against the reference. The paper describes three variant positions;
that discrepancy is unresolved. Main methods specify Rosetta E. coli and His-tag
removal with thrombin, while the RCSB summary names BL21(DE3). Exact assayed-batch
sequence identity remains a same-study approximation.

Chains A/B each contain observed author residues 36-373: 338 residues and 2,517
unique protein heavy atoms. Every observed standard residue has its expected
heavy-atom count; 18 declared N-terminal residues per chain are unmodeled. This
is a complete heavy-atom proxy for the observed region, not a claim that missing
residues were reconstructed.

The frozen selection rule is **assembly 1, lexical first protein author chain,
then lexical first attached glycan label chain**. This selects chain A's
N-terminal site, glycan G, independently of scores.

| Protein chain | Bound glycan label chain | Native Ca author ID |
|---|---|---|
| A | G **(primary)** | A:401 |
| A | H | A:402 |
| B | I | B:401 |
| B | J | B:402 |

The primary center is the arithmetic centroid of the 11 nonhydrogen atoms of
GLA G:2: **(21.664, 44.973364, 26.502545) Å**. The original depositor alias is
A/MLB501 in `_pdbx_branch_scheme`; the remediated `_atom_site` uses GLC G:1 and
GLA G:2. GLA has 11 atoms because its anomeric oxygen is replaced by the linkage
to GLC O6; only that O6 of the glucose tail is observed. The unresolved tail
provides no atomic pose reference.

Native Ca A:401, label chain S, coordinates GLA O3/O4 at 2.357/2.844 Å and
protein Asp97 OD1, Asp98 OD2, Asn138 OD1, Ser143 OG, and Pro144 O. The
[site audit](../data/archive/independent_specificity_sources/POL_site_audit.json)
records all four sites, exact distances, source hashes and heavy-atom checks.
Calcium/water handling, box extent, protonation and metal parameterization
belong in the computational freeze. This is a ligand-informed known-site
experiment; no blind-site claim follows.

## Retained alternatives and exclusions

| Candidate | Evidence recovered | Qualification decision |
|---|---|---|
| [Duff 2011, ConA](https://doi.org/10.1039/C0OB00425A) | New septanoside 7: Ka 840 ± 40 M⁻¹; 8 and 14: operational NB; explicit Table 1 buffer and concentration ranges | Promising; exact stereochemical drawings, C2 relationship and construct correspondence still require recovery. Pyranoside controls 9/10/11 are imported across earlier sources and are not promoted as one assay. |
| [Bouckaert 2004, FimH](https://doi.org/10.1111/j.1365-2958.2004.04415.x) | Mannose 2.3 µM, glucose 9.24 mM, galactose 100 mM in source tables | Quantitative weak-binding lead. Galactose has a finite table value, not an NB label. Exact competition doses and structure/construct correspondence remain open. |
| [Ficko-Blean 2009, NagHCBM32](https://pmc.ncbi.nlm.nih.gov/articles/PMC2937040/) | GlcNAc signal; tested GalNAc lacks a characteristic UV difference signal | Negative screen used an unspecified excess of solid sugar. Its concentration is not recovered. |
| [Weatherman 1996, legume lectins](https://kiesslinglab.com/sites/default/files/labs/kiessling/pdfs/1996.weatherman.biochem.pdf) | Actual finite and NB C-glycoside measurements | Negative monovalent 3/4 is an assayed diastereomer mixture; other relevant ligands are bivalent. No separately measured component labels are inferred. |

The archive contains six raw PDB/CCD files under the
[PDB archive CC0 policy](https://www.rcsb.org/pages/policies#usage), plus derived
geometry and a checksum manifest. Article bytes are not redistributed; their
public retrieval locations and verified hashes are recorded.
