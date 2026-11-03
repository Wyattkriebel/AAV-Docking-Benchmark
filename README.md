# AAV Docking Benchmark

A compact comparison of existing local AAV9 docking results: does using the
ligand-stripped holo structure or the independently solved apo structure change
the returned scores and agreement with the deposited modeled galactose?

**The six-run comparison is complete.** See the
[results and interpretation](comparison/RESULTS.md). Both states return favorable
Vina scores with substantial displacement from the chosen modeled reference.
These are descriptive results for one system; the coordinates have not been
qualified as atomic pose truth, and the comparison does not establish binding
specificity or general docking reliability.

## Run the comparison

Python 3.11 or later is sufficient; no packages, downloads or docking programs
are needed to analyze the included outputs.

```sh
python3 pilot/pilot-script.py
```

The command prints the comparison. To verify the committed results or regenerate
them from the archived inputs:

```sh
python3 pilot/pilot-script.py --check
python3 pilot/pilot-script.py --write
```

All six original pose files and logs are in [comparison/inputs](comparison/inputs).
[SOURCES.json](comparison/inputs/SOURCES.json) identifies the source snapshot,
checksums and the few extracted metadata fields. The analysis checks every input
hash, matches atoms by name, preserves returned pose rank and requires all three
configured seeds in each state. It reproduces the saved-output analysis;
preparation and native docking are not rerun by these commands.

## Scope

The two docking inputs each retain the same three whole subunits around a known
site. They are local neighborhoods, rather than full capsids. Search information,
software settings, receptor alignment and limits are recorded with the results.
The three seeds are paired computational repeats, not six independent systems.

The original [protocol](PROTOCOL.md) proposes a larger scientific study of site
discovery, pose recovery and assay-backed specificity. Those questions remain
open. This change supplies the smaller AAV9 comparison without requiring the
larger benchmark framework. Historical cropped-region [fpocket](pilot/fpocket)
and [P2Rank](pilot/p2rank) outputs are preserved as separate exploratory material;
their rankings are not pooled into this apo/holo docking comparison.

## Checks

```sh
python3 -m unittest discover -s tests
```

The small test suite checks atom-order independence, known geometric controls,
invalid pose rejection, complete run accounting and archive integrity.
