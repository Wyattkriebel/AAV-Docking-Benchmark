# Reproducing site localization and Gal9C specificity

There are two distinct operations: verify the retained evidence without native
predictors, or execute the frozen workflows again in separate directories.
A new execution is a replay of a development study whose outcomes are already
known. It does not become a held-out or newly prospective experiment.

## Verify the retained evidence

From the repository root, use the pinned Python environment:

```sh
python3.11 -m venv .venv
.venv/bin/python -m pip install -e '.[test]'
.venv/bin/python -m pytest -q
.venv/bin/python scripts/verify_site_predictions.py --output results/site_discovery/prediction_verification.json --check
.venv/bin/python scripts/verify_site_experiment.py --numeric-tolerance 1e-10 --surface-area-tolerance 1e-6
.venv/bin/python scripts/regenerate_site_report.py --check
.venv/bin/python scripts/regenerate_specificity_report.py --prepared-root results/gal9_specificity --check --numeric-tolerance 1e-10
```

The site prediction verifier restores the archived raw predictor files and both
frozen input revisions, checks their hashes, and reparses the original outputs.
The site experiment verifier recalculates the site endpoints from those inputs,
retained candidates and independently pinned native references. The specificity
report reads archived preparation files and all 36 raw docking outputs; it
verifies the commands, input fingerprints, atom maps, scores and log tables
before computing the assay comparisons. The explicit `1e-10` tolerance applies only to named derived geometry
values for cross-platform linear-algebra roundoff; identities, counts, labels,
decisions and archived bytes remain exact. Omit it to require exact numeric
recomputation. No fpocket, Java, P2Rank or Vina installation is required for these evidence checks. The full test suite includes
an optional live-fpocket test, skipped unless its executable is explicitly set.

Predictor input archives contain historical absolute paths as provenance. Their
coordinate, map and receipt bytes are verified after relocation; the verification
commands do not attempt to open those historical paths. Publication-source
checksums are declarations from the source audit, rather than a claim that all
publication bytes have been redistributed and reverified locally.

The default site replay requires exact values. The explicit cross-platform
flags above allow 1e-10 for named derived geometry values and, separately,
1e-6 Å² only for the accumulated probe-exposed surface-area total. Linux
showed a 4.05e-8 Å² accumulation difference in the approximately 859,673 Å²
capsid area. Coordinates, scores, ranks, counts, identities and success
decisions retain their stricter checks; the verification records the actual
per-field differences.

Specificity JSON checking can also explicitly allow 1e-10 only for the existing
derived-correlation and pKD fields. Raw energies, source measurements, labels,
counts, graphs and hashes remain exact; the Markdown report remains byte-exact.

## Execute the site workflows again

The retained fpocket build recipe and portable P2Rank runtime support **macOS
arm64**. fpocket needs the Apple command line compiler tools. P2Rank downloads
its pinned portable Java runtime, without changing the system Java installation.
A different native platform requires a reviewed build/runtime lock. The evidence
verification above remains usable on Linux and does not install these tools.

Run the following from the source checkout. Replace the two example replay paths
with new locations outside that checkout:

```sh
.venv/bin/python scripts/run_next_experiments.py install-site-tools
.venv/bin/python scripts/run_next_experiments.py site-prepare --work-dir /tmp/aav-site-replay --reprepare
.venv/bin/python scripts/run_next_experiments.py site-run --work-dir /tmp/aav-site-replay
.venv/bin/python scripts/run_next_experiments.py site-evaluate --work-dir /tmp/aav-site-replay
```

`site-prepare` restores both exact frozen receptor revisions and all native
source bytes from checksum-pinned local archives. `--reprepare` additionally
rebuilds every complete protein receptor directly from its source using the
current serialization revision, then requires exact equality of all author atom
identities and coordinates with the frozen inputs. These rebuilt files are
additional preparation controls; they do not replace the frozen method inputs.
The current serializer includes the terminal mmCIF comment delimiter that the
original fpocket reader required.

The recorded preparation smoke replay is preserved in
[preparation_replay.json](../results/site_discovery/preparation_replay.json).
It rechecks all nine actual receptor rebuilds against their frozen source/input
hashes and exact parsed author atom identities and coordinates, and all 12 Gal9
restorations against every declared archive file hash. Original and current
serialization bytes can differ; the receipt records both byte hashes and equal
atom-data hashes. It identifies the current command version and original replay
markers without claiming an execution-time script snapshot was frozen.

`site-run` executes the corrected fpocket input, the original P2Rank input, and
all 100 surface-baseline draws per condition. Each receptor is processed once
per predictor using the frozen flags. The surface exposure geometry is computed
once per receptor and reused across the frozen seeds. Capsid operators are
recomputed from receptor C-alpha atoms alone, checked against the previously
frozen source-only geometry, and used for symmetry-aware surface sampling.
Native coordinates enter only at the separate `site-evaluate` stage.

Each command writes solely inside its separate replay workspace. Preparation
requires a new or empty directory. Existing method output directories are
rejected, retaining earlier attempts and failures. To run methods in separate
steps, select each exactly once:

```sh
.venv/bin/python scripts/run_next_experiments.py site-run --work-dir /tmp/aav-site-replay --method fpocket
.venv/bin/python scripts/run_next_experiments.py site-run --work-dir /tmp/aav-site-replay --method p2rank
.venv/bin/python scripts/run_next_experiments.py site-run --work-dir /tmp/aav-site-replay --method random_surface
```

Tool paths may be supplied with `--fpocket`, `--p2rank`, and `--java-home`; the
adapters still enforce the locked installation and runtime. The default paths
point to the source checkout's verified local `tools/` installations. P2Rank's
per-receptor prediction timeout remains 600 seconds. fpocket retains its
3600-second operational timeout. Failures are recorded and execution continues
through the eligible cases; `site-evaluate` retains failures in its denominator.

The principal output is
`/tmp/aav-site-replay/results/site_discovery/evaluation.json`, with per-condition
reference and detailed evaluation files beside it. The historical Markdown
report is regenerated from the retained original study, using the verification
commands above. It includes the original nine serialization failures and the
observed corrected-run failure categories. Those historical observations must
not be copied into a new run as though they happened again.

fpocket's upstream source seeds randomness from wall-clock time and exposes no
command-line seed. Its native outputs and success/failure outcomes may differ
on replay; bitwise execution repeatability is not claimed. The archived raw
outputs, rather than a new execution, support reproduction of the reported
observations. P2Rank and the surface baseline use their frozen seeds.

## Execute Gal9C docking again

The Gal9C run uses the original 12 archived preparations, one shared receptor
and box, and three frozen seeds per compound. The box is informed by a prior
Gal9C glycan structure; this arm does not test blind localization or atomic-pose
recovery. [The configuration-field interpretation](../benchmarks/gal9_field_interpretation.json)
identifies inactive pose and eligibility fields inherited from the earlier
development template without changing any frozen bytes. Conditional source `NB` labels mean no binding detected at the tested
2 mM concentration and remain categorical, with no invented numerical affinity
bound.

```sh
.venv/bin/python scripts/run_next_experiments.py gal9-prepare --work-dir /tmp/aav-gal9-replay
.venv/bin/python scripts/run_next_experiments.py gal9-run --work-dir /tmp/aav-gal9-replay
.venv/bin/python scripts/run_next_experiments.py gal9-report --work-dir /tmp/aav-gal9-replay
```

If absent, the run command installs the frozen official Vina 1.2.7 executable
for a platform listed in `benchmarks/vina.lock.json`. A verified installation
can instead be supplied using `--vina`; the run/report validators still enforce
the allowed executable hash. Use the canonical executable filename
`vina-1.2.7`, which is part of the audited command convention.

The run executes 36 scheduled docking jobs with CPU 2, exhaustiveness 16, up to
20 returned modes, and the frozen 24 Å box. Existing case directories are
rejected. A failed run is retained, and the command returns failure after
processing all compounds. The report explicitly marks incomplete panels
ineligible for aggregate discrimination metrics. Successful panels produce
`/tmp/aav-gal9-replay/results/gal9_specificity/summary.json` and
`/tmp/aav-gal9-replay/docs/GAL9_SPECIFICITY_RESULTS.md`.

The replay copies the original plan and freeze receipt unchanged, alongside a
new `replay.json` identifying its historical development role. It never issues
a new pre-outcome freeze declaration. Native input-generation coordinates are
only atom-mapping checks; no RMSD against an independently generated input is
reported as an experimental pose endpoint.
