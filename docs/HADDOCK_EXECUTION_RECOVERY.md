# Preserving the first completed HADDOCK run and resuming orchestration

The first frozen main-cohort native job, `strict/2RDK/fpocket/box01/s104729`, generated all20 rigidbody and all20 flexref models normally. Its Python driver then failed while serializing the execution receipt: a documented relative CLI plan path was passed to `Path.relative_to()` against an absolute root. The frozen helper, plan, configurations and scientific inputs remain unchanged. No native rerun or new cohort condition is part of this recovery.

The original failed orchestration log remains `results/haddock_paired_cohort/execution_driver.log`. The first job's stdout/stderr, actual configs, PSFs, CNS logs, model files and I/O metadata remain in their original paths. `aav_benchmark/haddock_cohort_resume.py` is a separate wrapper. It normalizes the root, original plan and original freeze to absolute paths before calling the unchanged original `execute_job()` for future unique jobs.

## Recovered scientific outputs and missing telemetry

`orchestration_recovery/orphan_completion_proof.json` pins the preserved raw files and checks:

- Original frozen plan, helper, execution freeze, native source/runtime, config and all six input hashes match.
- The actual raw native configuration is byte-identical. Every explicit module parameter matches effective stage parameters. Both stage AIR copies are byte-identical.
- HADDOCK's input sanitization removes only serialization records such as CRYST1/TER: all760 protein and23 glycan atom identities, elements and coordinates are exact.
- Both actual PSFs parse; every one of40 generated models has the complete1083-atom inventory and finite coordinates.
- All20 initial and20 refined models have finite native scores and energy inventories. All20 parent relationships and inherited seeds match; refinement received the full initial pool.
- All40 docking CNS logs terminate normally. The native stdout has one normal final Finished marker and explicitly reports6minutes13seconds, or373seconds.

An independent read-only audit confirmed those facts. It also found the native stderr contained two NumPy empty-slice warnings, without a traceback. No experimental accuracy or self-analysis CAPRI result file was inspected. Native HADDOCK's automatic post-analysis uses its lowest-score model as a self-reference when no experimental reference is supplied; those raw files are retained, and are not treated as experimental pose accuracy.

The native subprocess return code was not captured in a retained receipt before the later Python exception. The recovered record therefore uses **`process_return_code:null`**, rather than inferring0 from completed products. Original process CPU and external measured wall time are also null. Its373-second `elapsed_seconds` is explicitly labeled as an integer-second native log report. Output hashes are anchored posthoc during recovery; they are not represented as contemporaneous subprocess telemetry.

After a separate coordinator recovery freeze, the wrapper writes a new `adoption_receipt.json`, then creates the previously absent `execution.json` with `status="complete_outputs_recovered"`, `native_rerun=false`, the complete original model inventory and a pinned `recovery_receipt`. The source-only completion proof and original failed driver log remain independently pinned. Existing terminal receipts are checked and never rewritten. An unreviewed native directory without a terminal receipt blocks resumption rather than being rerun.

## Future resources and bounded concurrency

There are38 remaining unique native jobs, with the original108-cell denominator unchanged. Every native job retains one CPU, the same20 initial model slots, inherited seeds, force field, AIRs and refinement settings. No failed model slot is replaced and no native job is submitted twice.

Each future job executes in its own spawned process worker. A separate `measured_resources.json` records `time.perf_counter()` wall time and `resource.getrusage()` user/system counters for the worker and its waited-on children before and after the call. This avoids attributing sibling jobs' CPU to one another. The resource scope and original OS counters are retained; worker startup before the measured call is excluded. Raw RSS counters are retained without treating their difference as a per-job memory peak. The [Python resource API](https://docs.python.org/3/library/resource.html) defines these process/child counters, and [ProcessPoolExecutor](https://docs.python.org/3/library/concurrent.futures.html#processpoolexecutor) provides process isolation. The wrapper uses one task per child process.

The coordinator's new freeze must declare the operational rule before resumption: at most two concurrent HADDOCK jobs while `results/blind_pose/execution_progress.json` has not finished all156 scheduled seed runs; up to six only when `scheduled_seed_runs=156`, `canonical_terminal_seed_runs=156`, `running_seed_runs=0`, and `status="complete"`. Missing or partially written status fails closed to at most two. Only workflow-completion fields are consulted, never accuracy. Exact status snapshots and effective concurrency changes are retained. This changes scheduling and measured wall cost, not model parameters or scientific candidate counts.

`concurrency_control.json` permits only `active_job_limit` and `stop_after_running`. Its requested limit is capped by both the new freeze and the old-workflow resource rule. Setting `stop_after_running=true` drains active jobs without submitting pending ones; a later invocation resumes from unchanged terminal receipts. No automatic model retry occurs. The native helper is no longer the overall coordinator.

## Review and execution boundary

The candidate recovery plan, completion proof and unauthorised freeze proposal are under `results/haddock_paired_cohort/orchestration_recovery/`. Preparation does not adopt metadata or launch another native job. Root must review and create a separate receipt with `status="frozen_before_HADDOCK_orchestration_resumption"`, `execution_authorized=true`, the exact recovery-plan SHA, CPU/concurrency cap and explicit resource rule. Only then may it run:

```sh
.venv/bin/python -m aav_benchmark.haddock_cohort_resume \
  --recovery-plan results/haddock_paired_cohort/orchestration_recovery/recovery_plan.json \
  --recovery-freeze PATH_TO_ROOT_RECOVERY_FREEZE.json \
  --control results/haddock_paired_cohort/orchestration_recovery/concurrency_control.json
```

Thirteen regression tests reproduce the original relative-path exception, prove absolute normalization, verify the frozen helper refuses the existing first job before a subprocess call, replay the actual complete products and telemetry gaps, reject missing models/energies/lineage/completion markers, test CPU counter scope, enforce the2→6 workflow-only resource rule, reject control changes to native CPUs, verify existing receipt immutability, require the separate recovery freeze, and prevent duplicate coordinator submission with an exclusive lock. No new native cohort job or experimental accuracy calculation is executed by these tests.

Portable reporting reads the retained receipts, raw models/PSFs/configs/AIRs/I/O and measured-resource sidecars. It does not invoke the native software validator or require the370 installed software files locally. Their historical metadata remains in the original public archive. Recovery concerns only the existing frozen main cohort.
