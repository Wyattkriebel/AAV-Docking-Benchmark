# Global report readiness gate

`aav_benchmark.global_pose_report_gate` adds a separately frozen coordinator around the existing report readers. It changes no frozen reader, native execution helper, input, score or metric definition. Its default operation returns workflow counters and computes no accuracy. No live probe or accuracy run was used to develop its synthetic tests.

All four gates must pass before the coordinator can call any unpatched accuracy reader:

| Gate | Fixed denominator | Completion requirement |
| --- | --- | --- |
| Main paired Vina/Vina-Carb | 36 conditions, 108 seed pairs, 216 method outcomes | All 216 outcomes terminal, including explicit failures |
| HADDOCK | 108 eligible cells; 39 native-ready cells | All 108 terminal, including source, assembly and unavailable-prediction failures |
| Assembly control | 12 control and 12 comparator conditions, 72 seed outcomes | All 72 outcomes terminal |
| Common-pose union | 108 pools, 216 source-method slots and 216 rescore-method slots | All source slots and score slots terminal; every retained pose has its required fixed score records |

These denominators are overlapping views, not independent samples. The assembly comparison includes main-cohort outcomes. The common-pose union reuses the main native searches. The older development-cohort report retains its own existing gate and is outside this four-way contract.

## Source-specific validation probes

Each probe runs in a separate Python subprocess. The three existing readers retain their full source, archive, execution and terminal-state validation. Immediately after each reader's pending branch and full-denominator check, its first reference-context call is replaced by a private readiness sentinel. The sentinel does not inspect its arguments. It finds the pinned entrypoint's caller frame, verifies the fixed denominator and full terminal counter in that frame, and raises an exception derived directly from `BaseException`.

The coordinator catches only that private sentinel. A recognized pending result is reduced to fixed denominators, its terminal count and the number of pending identifiers. A normal complete result, malformed counters, an unexpected exception, a timeout or entry into a blocked function fails closed. Pending results never release another reader's accuracy path.

Metric, fit, site, statistical and native-execution functions are blocked by name and by original function identity across loaded benchmark-module namespaces, including aliases. Statistic functions and process-launch APIs have additional identity blocking. Patches are restored on every exit. The fixed-union adapter supplies its own pure `probe_workflow` API; its probe additionally forbids entering reference context at all.

This is a guard for these reviewed, source-pinned readers, not a general sandbox for arbitrary Python. Independent source review must establish that the pinned reference-context call is the first metric boundary after each complete local gate. The AST record exposes that call location; it does not prove the reader's semantics. A changed reader requires a new reviewed contract rather than acceptance based on line numbers alone.

## Pins and authorization

The contract binds the running guard source, each actual metric plan and evaluation freeze, each reader's own source pin, its complete declared metric-source map, the reviewed context boundary and all fixed counters. Reader imports must resolve inside the reviewed repository. The union plan and receipt must both state that they were frozen before union accuracy. The existing three evaluation plans retain their original paths:

```text
benchmarks/structural_paired_method_evaluation_v2.json
results/structural_pose_v3/paired_evaluation_v2_freeze.json
benchmarks/haddock_pose_evaluation.json
results/haddock_paired_cohort/evaluation_freeze.json
benchmarks/structural_assembly_pose_evaluation.json
results/structural_assembly_pose_v3/evaluation_freeze.json
```

After independent review and the union evaluation freeze, create a proposal without invoking any input or metric reader:

```bash
.venv/bin/python -m aav_benchmark.global_pose_report_gate \
  --make-proposal benchmarks/global_pose_report_gate.proposal.json \
  --union-plan benchmarks/structural_paired_common_pose_evaluation.json \
  --union-freeze results/structural_pose_v3/common_pose_union/evaluation_freeze.json
```

Root must separately freeze the reviewed contract before live probes. That receipt has status `frozen_before_global_report_readiness_probes`, an exact `contract` path/SHA256 pin and `validation_only_probes_authorized: true`. Real generation additionally requires `real_reports_when_all_ready_authorized: true`. A source proposal is not authorization. The guard creates neither scientific freeze nor approval receipt itself.

The default CLI mode runs only the four blocked probes:

```bash
.venv/bin/python -m aav_benchmark.global_pose_report_gate \
  --contract benchmarks/global_pose_report_gate.json \
  --freeze results/global_pose_report_gate/freeze_receipt.json
```

## Explicit report generation

Only `--generate --output FRESH_DIRECTORY` can call the real metric readers. It first rechecks the source/plan pins and requires a fresh four-way readiness conjunction. It rechecks those immutable pins again, then calls the original three readers with pending disallowed. The union report receives the main paired report just generated in memory. No cached metric-file input is accepted by this coordinator.

The four reports are serialized with sorted keys and `allow_nan=False` before any output directory is published. They are written to a temporary sibling directory and renamed together to the requested fresh directory inside the repository. A failed calculation or serialization publishes no report directory. The original native outputs and frozen reader files remain unchanged.

The coordinator preserves method failures in the declared denominators. It neither retries native jobs nor repairs source failures. HADDOCK retains its guided AIR information distinction; readiness does not imply equivalence to a hard search box or to another method's input information.

## Synthetic validation

`tests/test_global_pose_report_gate.py` uses fake entrypoints, fake subprocess outputs and a tiny source-only pinned repository. It checks early-context rejection, alias blocking and restoration, exception propagation, fixed-count and source/plan tampering, all-four conjunction, pending isolation, fresh in-memory union handoff and strict atomic publication. No synthetic test launches a native job or invokes an actual accuracy reader.

```bash
.venv/bin/python -m pytest -q tests/test_global_pose_report_gate.py
```
