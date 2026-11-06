# Gal9C receptor sequence compatibility

Audit: 30 October 2026. **3NV3 passes the sequence compatibility gate for an
exploratory Gal9C CRD approximation.** Its complete 138-residue deposited
sequence is exactly human Gal9 canonical residues **218–355**, with no
substitutions. The reported assay construct spans 205–355, adding 13 N-terminal
residues. The entire annotated C-terminal galectin domain and every checked
binding-site residue are present in the crystal coordinates and identical.

This supports running the frozen Gal9C specificity cohort with a disclosed
construct approximation. It does not establish that the isolated crystal core
reproduces the full tagged assay construct. No ligand, assay label, or frozen
specificity record was changed during this audit.

The full sequences, residue-by-residue site checks, source checksums, and
decision fields are in
[`benchmarks/gal9_receptor_audit.json`](../benchmarks/gal9_receptor_audit.json).

## Exact sequence mapping

The [3NV3 coordinate deposition](https://files.rcsb.org/download/3NV3.cif)
contains one protein entity and chain A. Its entity sequence has 138 amino acids.
Direct comparison with the
[reviewed UniProt O00182 canonical sequence](https://rest.uniprot.org/uniprotkb/O00182.fasta)
gives an uninterrupted 138/138 match to residues 218–355. The 136 modeled
residues match canonical 220–355, also without substitutions.

| Sequence or coordinates | Inclusive range | Length |
|---|---|---:|
| Reported assay region, canonical O00182 | 205–355 | 151 |
| Declared 3NV3 entity, canonical O00182 | 218–355 | 138 |
| Declared 3NV3 entity, author numbering | 186–323 | 138 |
| Observed 3NV3 coordinates, canonical O00182 | 220–355 | 136 |
| Observed 3NV3 coordinates, author numbering | 188–323 | 136 |
| UniProt annotated C-terminal galectin domain | 227–355 | 129 |

The omitted assay N terminus is `GQMFSTPAIPPMM`, canonical 205–217. The first
two declared structural residues, canonical Y218/P219, are unmodeled. All 15
positions precede UniProt's annotated C-terminal domain. Every residue of that
129-residue domain is modeled. The first modeled residue, canonical H220, is
27.10 Å from the nearest native GAL heavy atom; this distance describes the
observed boundary and does not reconstruct the missing residues.

The current [UniProt record](https://rest.uniprot.org/uniprotkb/O00182.json)
is sequence version 2, last updated 15 July 1998. The assay's reported start is
also consistent with its published forward primer: the seven gene-specific
codons translate to `GQMFSTP`, canonical 205–211. The actual 2019 expression
batch or plasmid sequence was not recovered.

## Historical isoform names

The deposition calls its source a “short isoform variant,” references Q53FQ0,
and uses residues 186–323. Current Q53FQ0 requests redirect to the O00182
canonical record because Q53FQ0 is now a secondary accession.

An independent comparison to the current
[O00182-2 isoform sequence](https://rest.uniprot.org/uniprotkb/O00182-2.fasta)
also gives an exact match at 186–323. That isoform is currently named
**Medium** and deletes canonical 149–180, accounting for the 32-residue offset.
Current isoform 3 is named Short and deletes canonical 149–192. The deposition's
historical name therefore must not be equated to current isoform 3. Exact
sequence equality resolves the mapping:

`canonical residue = 3NV3 author residue + 32`

## Local binding-site check

The local audit uses protein heavy atoms within 6 Å of the independently
determined native GAL B3 heavy atoms in 3NV3. This selection uses no new docking
outcomes. All ten contacted residues are modeled and identical to the reported
assay region:

| 3NV3 author residue | Canonical residue | Closest GAL heavy-atom distance, Å |
|---|---|---:|
| R221 | R253 | 4.62 |
| H235 | H267 | 2.80 |
| N237 | N269 | 3.05 |
| R239 | R271 | 2.99 |
| V246 | V278 | 4.16 |
| R247 | R279 | 5.46 |
| N248 | N280 | 2.85 |
| W255 | W287 | 3.67 |
| G256 | G288 | 5.22 |
| E258 | E290 | 2.74 |

UniProt independently annotates canonical 267, 271, 281, and 287–293 as
C-terminal galactoside-binding positions. All ten positions are observed and
identical, including those outside the chosen 6 Å local shell. There are no
substitutions anywhere in the shared 138-residue sequence, so this result does
not depend on selecting one contact cutoff.

## Fusion and construct limits

The [2010 primary structure paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC2978625/)
explicitly describes expression as a GST fusion followed by overnight thrombin
cleavage and elution of the cleaved Gal9C target. The structural GST fusion was
therefore removed.

The assay's inherited
[Delaine et al. method](https://www.pure.ed.ac.uk/ws/portalfiles/portal/27411623/Delaine_et_al_2016_ChemBioChem.pdf)
reports Gal9C residues 205–355 encoded with thioredoxin in pET-32 Ek/LIC, from
I.M.A.G.E. clone 2208156. The recovered method does not state whether thioredoxin
was cleaved. Mahanti's 2019 SI adopts these methods without repeating an
assayed batch sequence or tag record. Retain those limits, the 13-residue
terminal truncation, and the two unmodeled residues in the run interpretation.

The sequence audit finds **no binding-site identity mismatch requiring an
unsupported-target failure**. It supports the known-site CRD comparison while
preserving the full-construct approximation. It does not qualify an atomic
ligand pose, establish a blind binding-site test, or supply a new receptor
family held out from earlier Gal3 work.
