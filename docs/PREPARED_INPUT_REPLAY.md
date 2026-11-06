# Replay the fixed preparation inputs

This bundle publishes the frozen preparation library and the immutable inputs
for six full-assembly cases and two source-chain A controls. It packages the
existing input bytes after execution began; it does not change the original
freezes or select inputs using docking outcomes. The preparation indexes retain
all 48 conditions: 28 prepared boxes, 18 preparation failures and two conditions
with no predicted candidate.

With the pinned repository dependencies installed, run:

```sh
python scripts/verify_prepared_inputs.py
python -m pytest -q tests/test_prepared_input_bundle.py tests/test_heldout_pose_inputs.py tests/test_heldout_assembly_profiles.py tests/test_heldout_assembly_pose_inputs.py
```

These commands require no `tools/`, `data/processed/`, network access or native
executables. The six actual ligand export cases and both monomer source profiles
are mandatory archive fixtures. Missing fixtures fail the tests instead of
skipping their replay. The export tests restore the original V2 chemical map
and V3 independent SDF/report into temporary directories; the existing conformer
is copied without a new embedding or minimization.

`results/prepared_input_foundation/prepared_inputs_manifest.json` pins a
deterministic archive of 698 existing files: the two preparation indexes, all
48 condition records and their prepared/failure files, and the 16 completed
predictor input streams that supplied their original ranks and centers. The
predictor files are retained pose-preparation prerequisites. No site or pose
accuracy is calculated from them here. Native docking outputs and live progress
records are excluded.

The verifier checks the packaging manifest, all archive/member hashes, sizes
and exact member identities, eight frozen producer hashes, four source-profile
archive pairs and the listed source/freeze documents. Archives are read in
memory; selected source CIFs are copied to fresh temporary files for exact
neighborhood reconstruction. Historical absolute paths resolve only against the declared
original repository prefix and then used as archive keys; those original paths
need not exist. The immutable scientific configurations are preserved byte for
byte.

The source checks require identical full-parent ligand SDF, PDBQT, canonical
map, chemistry and initializer report bytes across paired profiles. Repair
retains the original observed protein map. Monomer controls retain the full
observed assembly and map alongside the exact chain A subset. Modeled atoms
remain separately labeled source additions.

For all 28 prepared boxes, replay checks the original candidate rank/center,
the full canonical ligand identities, exact rigid-translation PDBQT bytes,
and exact whole-chain protein neighborhood PDB bytes under the frozen 8 Å
rule. It also checks that receptor PDBQT preparation preserved every input
heavy-atom identity and coordinate. These are input-correspondence checks;
there is no ligand fit, experimental recovery threshold or method comparison.
All 20 failed or empty input conditions remain represented by their original
condition and attempt records.

The preparation APIs and their historical CLIs are included with their frozen
bytes. Commands that generate atoms, run predictors or execute docking still
require their explicit scientific freeze, separately installed pinned software
and fresh output directories. The portable replay command does not invoke
those operations or reverify installed native software. Archived software and
reader receipts retain their historical provenance.

Successful output contains counts of archived inputs and retained conditions.
This input publication makes no docking-performance, binding, population or
model-training novelty claim. Scientific reports and native outcome archives
are published separately after their declared completion gates.
