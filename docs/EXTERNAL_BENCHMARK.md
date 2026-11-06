# Released examples and the exploratory Gal3 panel

The external cases validate implementation behavior. They do not constitute a
held-out or independently density-qualified accuracy cohort. The source
[HADDOCK protein–glycan dataset](https://github.com/haddocking/protein-glycans/tree/344e4a41a300b68b7ee3585895a460c863545ddb)
is pinned at commit `344e4a41a300b68b7ee3585895a460c863545ddb`; the
[associated paper](https://doi.org/10.1021/acs.jcim.4c01372) evaluates different
methods and information regimes. Reusing its prepared example with Vina is
not a reproduction of its published HADDOCK success rates.

## Frozen development inventory

`benchmarks/external_cohort.json` records the complete input inventory before
these new docking runs. It includes exact URLs and SHA256 digests for released
bound receptors/ligands, CCD monomers and the upstream Apache-2.0 license.
Their original bytes are also preserved in deterministic gzip archives under
`data/archive/external/` (about 0.2 MB total). Preparation restores and verifies
these archives before considering a network download; upstream revisions
cannot silently replace the frozen sources.
`benchmarks/sources/glycan_connectivity.json` separately records deposited
branch links and author/label identifiers. No bonds are inferred from ligand
atom proximity.

| Case | Input chemical object | Development limitation |
|---|---|---|
| 154L | NAG trisaccharide, 43 heavy atoms | Source reproduction; ligand density and validation concerns require independent review |
| 1GNY | Xylopentaose, 46 heavy atoms | Explicitly reviewed 40.35 Å longest box side; preserves 8 Å margins |
| 1KJL | LacNAc, 26 heavy atoms | Ligand-stripped isolated galectin-3 carbohydrate-recognition domain |
| 1QFO | Nominal GLC/GAL/SIA trisaccharide | Released reducing-end C1 handedness disagrees with CCD GLC; preparation fails visibly |
| 2J1T | Lewis-Y tetrasaccharide, 46 heavy atoms | Released bound receptor; biological assembly/local density not independently certified |
| 1PWB | Maltose, 23 heavy atoms | Native calcium retained with the fixed ion model described below |
| 2VUZ | Resolved seven-residue fragment, 90 heavy atoms | Not a full atomic reference for the nonasaccharide named in the entry title; native calcium and one unresolved protein residue |

The current parser verifies an exact heavy-atom set against CCD chemistry after
only the declared glycosidic leaving groups are removed. It checks each CCD
tetrahedral center against named-neighbor handedness in the deposited pose.
This catches the 1QFO anomer inconsistency without silently changing a chemical
identity. All eligible source cases remain in the preparation denominator.
An unsupported input is not removed to improve the conditional docking rate.

For the structural diagnostic, bound coordinates identify chemical
stereochemistry and the oracle box. They are discarded before ETKDGv3
embedding (seed 20261024, one thread) and converged MMFF94 optimization. The
fixed second embedding attempt uses random initial coordinates, the same
seed and a 1000-iteration cap; both attempts are logged when the first fails.
This development fallback was adopted after the initial 2VUZ embedding
failure, before any 2VUZ docking; it is not a claim of a prospectively frozen
confirmatory preparation procedure.
No native ligand coordinates or torsions seed the docking conformer. Meeko
heavy-atom identities are checked after PDBQT conversion.

1GNY's frozen x-span is 24.35 Å. The 8 Å margin on each side gives a 40.35 Å
box. This 0.35 Å exception to the automatic 40 Å review gate was reviewed
before that case's docking. It is an explicit development decision, recorded
in the machine-readable case definition.

`prepare_development(raw_dir, output_root)` in `aav_benchmark.external` writes
a case-by-case preparation manifest, including failures. Individual commands
use `fetch_external` and `prepare`; docking and scoring use the repository's
shared runner. Failed preparation attempts invalidate any earlier success
manifest. `chemistry.json` preserves successful embedding and atom-map checks
even if receptor preparation fails later.

## Native calcium and unresolved sequence gaps

1PWB and 2VUZ initially failed because the released receptors use `CA2:CA+2`
ion names. The development extension retains these original identities and
coordinates, parameterizes an isolated CCD-verified `[Ca+2]` through Meeko,
and appends its single `Ca` atom with explicit charge +2 to the rigid receptor.
The full original receptor heavy-atom set and coordinates are checked again.
This is a fixed Vina ion model; it does not model specialized coordination or
polarization and has not been validated against a calcium-dependent assay.

2VUZ's prepared source also renumbers a genuine unresolved residue away:
source GLY90 and ASN91 correspond to deposited author GLY90 and ASN92.
`benchmarks/sources/2VUZ_receptor_numbering.json` records the exact atom-name/
coordinate mapping for all 128 protein residues, with the deposited CIF source
and its checksum. Restoring author numbering for Meeko preparation exposes
the gap at author residue 91. No loop atoms or peptide bond across the gap
are created. All 990 original receptor heavy atoms, including calcium, remain
at their input coordinates. The output PDBQT restores the original source
identifiers so contact evaluation uses the same atom identities. The preparation
log declares both identity frames and the unresolved sequence gap. This case
remains a partial-structure development diagnostic.

These extensions were made during development, before docking either newly
supported case. `benchmarks/development_preparation.json` retains the initial
blocked preparation history beside the latest supported-case counts.

## Assay-qualified quantitative ranking diagnostic

`benchmarks/galectin3_itc_panel.json` captures three measured solution ITC
affinities for human galectin-3 from
[Sindrewicz et al. (2019)](https://doi.org/10.1038/s41598-019-47658-8): LacNAc,
lactose and galactose. The ledger retains the reported estimates, errors,
construct, temperature, buffer and concentration information, and marks
unreported conditions as unknown. All three compounds bind. The panel cannot
measure binary specificity or a false-positive rate.

The assay uses full-length Gal3, whereas the 1KJL docking receptor is an
isolated domain. This is a disclosed domain approximation. Every ligand uses
the same LacNAc-informed 1KJL box and receptor. The reducing end of each free
sugar can mutarotate: alpha and beta forms are built separately from their
distinct CCD identities (NDG/NAG, GLC/BGC and GLA/GAL). Their weights are not
assumed or fitted. The nonreducing galactose linkage remains beta.

The analysis was frozen before scoring in
`external_cohort.json:gal3_panel_frozen_analysis`: within each anomer, average
the seed-specific Top-1 Vina energies, compare the three ligand ranks with
measured affinity ranks, and show a heavy-atom-count baseline. Report
Spearman/Kendall descriptively. Three ligands and one target do not support
held-out generalization or a family bootstrap. The LacNAc/lactose comparison
has equal charge and differs by three heavy atoms; it is not an exact size
match.

`prepare_gal3_panel(raw_dir, receptor_dir, output_root)` writes six prepared
directories. Their `reference.json` contains independent input coordinates
solely for atom-identity verification and has `pose_reference_eligible=false`.
Those coordinates are not experimental truth and must never be used for
native-pose RMSD. No binding free energy or KD is inferred from Vina scores.

## Offline controls and source attribution

Run `.venv/bin/python -m unittest discover -s tests -p test_external.py -v`.
The controls use a small byte-identical subset of pinned source ligands and
CCD definitions in `tests/fixtures/external/`. They test complete atom maps,
deliberate atom loss, opposite anomers, preservation of linkage stereochemistry,
reference-independent embedding, duplicate/alternate atom rejection and
source checksum refusal. The fixture's `HADDOCK_DATA_LICENSE` is the original
unmodified Apache-2.0 license; source URLs and checksums are in the cohort
manifest. These tests validate chemistry handling, not biological accuracy.
