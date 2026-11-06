# Reference-independent site-discovery methods

This is a development extension of `PROTOCOL.md`. A site-discovery run consumes
the complete supplied protein receptor and frozen method settings. Native
ligands, experimental contact lists, oracle boxes and reference-derived
neighborhoods are absent from discovery inputs and method invocations. Native
coordinates are used only by the separate evaluation stage after predictions
are retained. Successful execution does not qualify an external receptor's
biological assembly or establish held-out performance.

## Receptor input

`aav_benchmark.sites.prepare_site_receptor` accepts one source model and writes a
protein-only heavy-atom receptor. It retains every source protein chain and
canonical amino-acid heavy atom, including chains distant from any known site.
There is no radius, center, ligand identifier or chain-selection argument.
Hydrogens, water, ions, metals, glycans and other nonprotein residues are
removed uniformly by chemical category. Alternate conformers, nonpositive
occupancy, duplicate identities, noncanonical/HETATM amino acids and nonfinite
coordinates require explicit curation and fail preparation. No protonation,
repair, minimization or alignment occurs.

Whole water and metal residues are classified before the PDB record type;
`ATOM`-encoded `CA2` ions in the external sources are still removed and counted
as metals. A protein alpha-carbon named `CA` remains a protein carbon atom.

The metal-free policy is shared by all site methods in this development arm.
It differs from the prior docking preparation that retained native calcium;
metal-dependent site recognition therefore remains a limitation. In particular,
the upstream fpocket heteroatom allowlist is not a general retained-metal
policy: calcium can be omitted internally. Removing every metal explicitly
avoids falsely claiming a retained-calcium experiment. This is not evidence
that ions are unimportant, and a retained-metal sensitivity arm would require
a separately frozen, audited implementation.

For 7MTZ and 7MT0, `assembly_id="1"` invokes the existing reviewed assembly
audit: one identity materialization and 60 distinct protein subunits. It never
uses the previous ligand-derived neighborhood. Other input receptors retain
their entire deposited model; the manifest explicitly does not certify their
biological completeness. The third-party benchmark receptor selection and
its original curation provenance remain qualifications on those cases.

Author chain, residue number, insertion code, residue name and atom name are
the stable atom identity. A deterministic compressed TSV also records original
label chain/sequence identifiers, serial, element and coordinates. The
preparation manifest hashes the source, final receptor, map and method lock,
records omissions and verifies the final heavy-atom identity set and maximum
coordinate displacement (at most 0.001 Å for serialization). Generated mmCIF
uses numeric entity identifiers for BioJava compatibility; where a PDB source
has no label sequence identifiers, these are assigned sequentially. Author
identifiers remain unchanged and missing original label IDs remain explicit
in the source map.

Use mmCIF for the full AAV capsid. 7MTZ has 247,860 protein heavy atoms, beyond
the decimal atom-serial capacity of legacy PDB. The upstream fpocket PDB reader
uses decimal `atoi`, so hybrid-36 output would lose serial meaning. The
preparer rejects PDB serial, chain, residue and coordinate overflow; it never
truncates an assembly to fit PDB. The shared mmCIF receptor can also be read
directly for the solvent-accessible surface baseline.

### Serialization compatibility revision 2

The initial nine scientific fpocket invocations failed immediately at input
reading, while P2Rank accepted the same valid Gemmi-produced mmCIF files.
These attempts are retained as implementation failures, not biological
non-discovery observations. Gemmi's serialization ended the atom loop at EOF.
The pinned fpocket binary's bundled VMD reader instead requires a line
beginning with `#` to end its atom-row scan. Its published
[parser source](https://www.ks.uiuc.edu/Research/vmd/plugins/doxygen/pdbxplugin_8C-source.html)
shows that EOF before this comment delimiter returns zero atoms. Appending the
delimiter resolves the observed failure without modifying atom fields,
algorithm source, method parameters or native-site information.

Revision-2 preparation therefore writes an explicit terminal `#` comment.
`revise_site_receptor_serialization(old_cif, new_cif)` handles already frozen
inputs: it requires a fresh destination, retains the original bytes and
manifest, copies the original atom map unchanged and verifies exact equality
of all parsed protein coordinates and author atom identities. A separate
revision receipt records both receptor byte hashes, the prior preparation
manifest hash and the unchanged atom-map hash. New scientific attempts must
follow a newly frozen input-revision receipt; earlier attempts and their costs
remain in the execution ledger. The runner copies the revision receipt into
the attempt directory and records its input serialization revision.

A paired software compatibility probe is retained in
`results/site_methods/software/fpocket/cif_compatibility_probe/`. Its independent
100-atom synthetic receptor failed through the actual adapter without the
delimiter and completed after adding it, with all 100 input atom identities
and coordinates recovered in fpocket's combined output. These artificial
pockets are excluded from every scientific denominator and site result.
The opt-in `test_real_fpocket_adapter_reads_100_atom_prepared_cif` also invokes
the real pinned binary on freshly prepared revision-2 input; enable it with
`AAV_FPOCKET_EXECUTABLE` set to that installation's executable path. A
serialization-only unit test verifies that revision 1 cannot be overwritten
and that the source atom map and coordinates remain exact.

## fpocket

The implementation uses official [fpocket release 4.2.3](https://github.com/Discngine/fpocket/releases/tag/4.2.3),
commit `4bb0d8447f62fee77e2c3c29f54b5fcaf5e2c066`. The exact codeload source
archive has SHA256
`62807fc21854e26a39e440e499cb27e1128a614509513921b83c3e239c35ae54`.
The root MIT license and bundled Qhull notice are checksummed in
`benchmarks/site_methods.lock.json` and copied alongside the build evidence.
The source archive also pins the shipped molfile plugin, which carries
University of Illinois copyright notices. No binary distribution is added to
the repository.

`install_fpocket` verifies that archive, extracts source/headers and the native
Apple arm64 plugin, and compiles `bin/fpocket`. The build recipe removes the
unsupported Apple arm64 profiling flag `-pg` and supplies the Qhull include
directory for that direct build target. Algorithm sources are unmodified.
The build manifest retains the compiler, exact build command, platform,
executable hash, stdout/stderr and raw help output. This builder is explicitly
limited to macOS arm64; another platform needs a reviewed build recipe.
The release source's help banner still says `fpocket 4.0`. The release,
commit/archive hash and compiled executable hash identify this installation;
the inconsistent banner is retained rather than relabeled.

These exact arguments are frozen before the scientific run:

```text
fpocket -f receptor.cif -m 3.4 -M 6.2 -i 15 -D 2.4 -A 3 -p 0.0 -v 300 -C s -e e -w m
```

They make the defaults in the pinned
[parameter header](https://github.com/Discngine/fpocket/blob/4bb0d8447f62fee77e2c3c29f54b5fcaf5e2c066/headers/fparams.h)
explicit: alpha-sphere radii 3.4–6.2 Å, at least 15 spheres, clustering distance
2.4 Å, three apolar neighbors, minimum apolar proportion zero, 300 volume
Monte Carlo iterations, single-linkage clustering with Euclidean distance and
mmCIF atom output. Ligand/pocket definitions, chain filters and energy grids
are not supplied. The settings were not selected from native-site recovery.

The unmodified upstream algorithm uses `srand(time(NULL))` for Voronoi
perturbation and Monte Carlo volume and provides no command-line seed. One
scientific execution per system is frozen. Start/end times and all raw
outputs are retained; bitwise repeated-run equivalence is not claimed. An
excluded installation/parser probe predates the freeze and is not included
in scientific denominators or evaluation. Re-running a failure creates a new
execution directory and must retain the first attempt and its cost.

Candidate rank preserves the upstream pocket number, which follows fpocket's
original pocket-score order; there is no native-distance or druggability
reranking. Every raw candidate is retained for later symmetry-aware merging.
The candidate center is the unweighted arithmetic mean of all alpha-sphere
centers in `pocketN_vert.pqr`, using the output's three-decimal Å coordinates.
This matches the alpha-sphere barycenter endpoint described in the pinned
[upstream guide](https://github.com/Discngine/fpocket/blob/4bb0d8447f62fee77e2c3c29f54b5fcaf5e2c066/doc/GETTINGSTARTED.md).
The footprint contains the receptor atoms contacted by those spheres, from
`pocketN_atm.cif`; a residue is included when it has at least one footprint
atom. The parser verifies every footprint atom's original author identity and
coordinate against the complete method input. It validates sphere counts and
rank continuity. The combined receptor output is also audited against the
entire input atom set, so an upstream-dropped chain or protein atom is a failed
execution. When fpocket reports no pockets it writes no combined receptor,
and this additional internal-reader audit is explicitly unavailable.
fpocket's atom files contain unquoted descriptor prose before
their atom loop; the adapter preserves the raw files and parses only the
`_atom_site` loop rather than treating that prose as valid mmCIF.

`run_fpocket` requires the preparation sidecar and pinned installation manifest
and rejects mutated inputs or binaries. It copies the exact prepared receptor
into a fresh run directory, records the command and hashes, and retains raw
stdout, stderr, descriptor/vertex/atom files and normalized predictions. Its
execution record includes elapsed time, child CPU usage and the process-family
child-memory high-water mark; the latter is explicitly not isolated per-run
RSS. An exit code of zero without output is a failure, because fpocket can
return zero after input errors. Its explicit `no pockets found` result is a
completed zero-candidate observation. Failed executions remain in the ledger
with their logs and error category.

The normalized candidate schema provides `rank`, `center_A`, `score`,
`alpha_sphere_count`, `residue_keys` (`chain`, number, insertion code, residue
name), `atom_ids` (the same plus atom name), and raw artifact references.
`center` and structured `residues` are compatibility aliases.

## Use and validation

```python
from aav_benchmark.sites import install_fpocket, prepare_site_receptor, run_fpocket

installation = install_fpocket("tools", evidence_dir="results/site_methods/software/fpocket")
prepare_site_receptor("data/raw/7MTZ.cif", "data/processed/sites/7MTZ/receptor.cif", assembly_id="1")
result = run_fpocket("data/processed/sites/7MTZ/receptor.cif", "results/site_methods/7MTZ/fpocket",
                     executable=installation["executable"])
```

The cohort execution/evaluation layer must freeze its case manifest and checksums
before invoking these scientific runs. Run failures, zero candidates and
reference-qualification failures are different outcomes and keep their
eligible denominators. Discovery does not imply successful pose docking.

`tests/test_sites.py` verifies ligand-coordinate invariance, whole-receptor
retention, uniform metal/water removal, source atom maps, numeric entity IDs,
author insertion codes, PDB overflow rejection, alternate/duplicate rejection,
the alpha-sphere barycenter, footprint identity/coordinate audits, mutated-input
rejection, nonzero-exit evidence, false-success rejection and explicit
zero-candidate handling, including ATOM-encoded calcium and full-output atom
loss. Synthetic output fixtures and subprocess doubles are
test controls only and are never scientific predictions.
