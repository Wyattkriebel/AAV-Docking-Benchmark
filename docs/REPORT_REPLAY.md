# Comparing regenerated report JSON

`scripts/compare_report_replay.py` compares two existing files. It does not run
search, scoring, reference qualification or accuracy metrics, and does not
replace the global workflow gate or input/source checksum verification.

The default requires identical bytes after validating both files as UTF-8 JSON:

```sh
python scripts/compare_report_replay.py saved/main_paired.json replay/main_paired.json
```

For a separate, explicit numerical portability check:

```sh
python scripts/compare_report_replay.py saved/main_paired.json replay/main_paired.json --numeric-tolerance 1e-10
```

The tolerance is absolute, with no relative term, and must be finite and between
0 and `1e-8`. The existing archive verifier uses the same maximum bound. This is
an implementation comparison option, not a change to any metric, success
threshold, scientific decision or uncertainty interval.

The field policy, `global_report_derived_paths_v1`, is defined by complete paths
in [report_replay.py](../aav_benchmark/report_replay.py). It permits roundoff only
for these reviewed output sections of the main, HADDOCK, assembly and union
reports:

| Output section | Derived values eligible for explicit tolerance |
|---|---|
| Retained poses and copied union pose metrics | Receptor-frame RMSDs, minimum protein distances, contact/interface fractions, bond/ring distance changes and linkage distances |
| Pose-rank summaries | Top-1, Top-5 and best-retained RMSDs |
| Protein-fit controls | Derived global/local backbone RMSDs, rotation/translation components and the computed translation-control RMSD |
| Site evaluation | Center-to-native distances, qualified-copy distances and coverage fractions |
| Conditional statistics | Recovery means, paired differences and interval endpoints |

Matching a field name elsewhere does not authorize a tolerance. Native method
and ranking energies, resource telemetry, raw coordinates, input centers,
confidence levels, thresholds, hashes, identifiers, integer counts, atom maps,
boolean decisions and every unknown path remain exact. Analytic zero native
RMSD and the known translation vector also remain exact.

Structural mode interprets JSON floating values as finite Python binary64
values. It ignores equivalent number formatting, whitespace and object-key
order; use the default byte mode when those bytes must also match. JSON integer,
floating, boolean, string and null types remain distinct. Signed zero remains
exact outside an explicitly tolerated derived path. Duplicate object keys,
nonfinite numbers and malformed or excessively nested input are rejected even
when the two files have identical bytes.

Tolerance bounds are checked using exact ratios of the represented binary
values. This avoids rounded subtraction admitting a difference just above the
requested bound. Diagnostics include the path, maximum difference, exact
difference/tolerance ratios, number of tolerated differences, file hashes and
sizes. `--differences-limit` controls only how many example differences are
printed; every value is still checked.

Lists are compared position by position without sorting. Pose/family identities,
ranks and integer atom maps therefore protect row order. A permutation of
derived vector values whose differences are already below tolerance is
indistinguishable from numerical roundoff; the comparator does not claim to
infer such a permutation. Raw coordinate vectors have no tolerance.

`--structural` requests exact parsed types and values while ignoring formatting,
without enabling numeric tolerance. A successful comparison exits zero;
rejected comparisons exit one and file-reading errors exit two. The module and
CLI require only Python's standard library and work with `python -S`.
