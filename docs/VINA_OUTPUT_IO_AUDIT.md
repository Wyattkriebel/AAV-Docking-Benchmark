# Auditing Vina's logged modes and saved pose prefix

AutoDock Vina 1.2.7 can log 20 modes and save fewer models. Its writer stops when a mode's energy exceeds the best energy plus `energy_range`. A valid saved file therefore contains the corresponding complete rank prefix. Requiring equal logged and saved counts rejects some valid output files.

The standalone `aav_benchmark.vina_output_audit` module checks this I/O behavior. It uses the Python standard library and reads only execution metadata, MODEL ranks, score tokens and complete-file hashes. It does not parse atom coordinates, open preparation or experimental-reference files, calculate accuracy, launch native software, change eligibility, or modify the existing frozen readers.

Run the read-only audit from any working directory:

```sh
python /PATH/TO/REPO/scripts/audit_vina_outputs.py --root /PATH/TO/REPO
```

The command prints JSON to stdout and returns 1 if a raw output or manifest fails validation, or if there are no successful outputs to audit. Running source records remain explicit pending records. Repeat `--scope` to restrict the audit, for example `--scope development --scope pol_specificity`. The default fixed allowlist covers the existing main experiments, controls and clean-environment replays; historical recovery copies, Vina-Carb and supplemental assembly outputs are excluded. Symlinked scopes, manifests and raw artifacts are rejected before opening them. Counts are run records, including replay copies, rather than independent experiments.

Each successful record includes its raw manifest, log and saved-file SHA-256 hashes and byte counts. Validation requires the Vina 1.2.7 banner, ordinary verbosity 1, an explicit 10 kcal/mol window, ordered contiguous log ranks, a complete saved rank prefix, and compatible score tokens. It checks one shared best energy across the whole prefix and suffix.

Ordinary Vina logs use four significant digits; PDBQT scores use three fixed decimal places. The validator calculates decimal rounding intervals, intersects log and saved intervals, and tracks nearest-even endpoint inclusion. Power-of-ten intervals are asymmetric: logged `10` has bounds 9.9995 to 10.005. Logged zero is exact under significant-digit formatting. An omitted suffix must be compatible with the writer's strict `energy > best + 10` rule. The report distinguishes a **proven outside-window** omission from a **boundary-ambiguous, writer-compatible** omission; ambiguity is never counted as proof. Unsupported verbose formatting, malformed tables, missing ranks, incompatible scores, and proven inside-window omissions fail validation. These conventions follow the [Vina 1.2.7 search log](https://github.com/ccsb-scripps/AutoDock-Vina/blob/v1.2.7/src/lib/vina.cpp#L493-L542) and [pose writer](https://github.com/ccsb-scripps/AutoDock-Vina/blob/v1.2.7/src/lib/vina.cpp#L915-L933).

The intervals describe continuous decimal rounding compatibility. They do not prove that every compatible value is an exact binary-double value, or authenticate a native execution. A finite exponent envelope bounds parsing; displayed log tokens can exceed the underlying double's range after rounding. Raw artifact hashes and recorded execution provenance remain separate evidence.

`tests/fixtures/vina_output_io_controls.json` contains synthetic text controls with no atom or coordinate records. `results/vina_output_io_audit/software_validation.json` records their hashes, replay results and exact software-file pins. This artifact is independent of molecular sources and native docking software. The tests also replay three already archived 20-logged/19-saved controls; no new scientific outputs are created.

In a clean checkout of the PR7 base plus the six files in this change:

```sh
python -m pytest -q tests/test_vina_output_audit.py
python -S scripts/audit_vina_outputs.py --scope development --scope pol_specificity
```

The second command demonstrates that the audit CLI needs no site packages or native executables. The first command requires pytest; the module itself has no third-party dependencies. Neither command changes source data. The current main experiments can still have equal logged and saved counts, so their successful I/O checks do not exercise the omitted-suffix case; the archived and synthetic controls supply that coverage.
