# Paired execution work sharing: proposal only

No live coordinator, reservation, process signal or native job was created by
this work. The implementation supplied here is a synthetic state model and
tests. It is **not ready for live execution**. The preserving coordinator and
its real process/log/output validators must be implemented, tested and reviewed
before a separate root operational freeze and before the first reservation.

Waiting until all 156 old jobs are independently verified terminal before
adding two CPU-2 workers is a reasonable resource gate. It would permit at most
three simultaneous paired searches while the original CPU-2 driver remains
active, and at most two after that driver exits. It does not change the fixed
36 conditions, three seeds, two methods, 216 eligible method-seed slots,
114 potential full searches before CHI failures, scientific parameters or input
bytes. No condition or receptor profile may reuse another profile's execution.
The protected supplemental arm is outside this proposal.

## Prefer a post-reservation log barrier over process suspension

The frozen driver
`results/structural_pose_v3/paired_execution_driver.py`, SHA256
`1b846c076983776c1a6aa251aed709a5b79f6283ff68c74383f59cb6dbd0d36c`,
prints and flushes `START_MODERN condition_id` before calling the original
modern engine. That engine refuses a fresh run when `runs.json` already exists,
before launching a search. Its native version query precedes this guard; it is
not a search.

The narrow no-signal protocol is:

1. Independently verify the old156 completion gate, root operational freeze,
   complete source/input bindings, exact legacy PID/start/command/source identity
   and direct regular stdout-log binding. Open that same log file descriptor for
   both frontier reads. Acquire one of two worker slots for the entire condition.
2. Read the latest complete `START_MODERN` prefix. Choose a strictly future,
   prepared condition in reverse fixed-manifest order. Never claim the current
   or an earlier condition. Never claim a preparation/predictor failure branch:
   those branches bypass the modern engine's manifest guard.
3. Through the one coordinator, check all modern and three VC destinations.
   Skip a fully verified terminal condition. Stop on an unowned empty manifest,
   partial output, unknown writer or unverified existing attempt. Do not retry
   any completed failure. Record an exclusive ownership intent, then create
   `runs.json` containing exactly `[]` through `O_EXCL`; never replace a file.
4. **After exclusive creation, reread the same bound stdout descriptor.** If
   the parent has reached this condition, abort without any worker native call
   and retain the reservation/receipt. If the log changed identity, truncated,
   restarted or regressed, stop. A first read alone is insufficient.
5. Only if the second frontier is still earlier, recheck the owned empty inode,
   bytes and receipt, all absent VC destinations and the source bindings. The
   same worker may then call the unchanged modern engine with `resume=True`
   on that owned empty manifest, followed by the unchanged three VC seed jobs.
   Hold the worker slot throughout. Record every actual result and cost under
   the same frozen condition/seed contract, plus separate operational receipts.

The exact duplicate-prevention invariant is:

> Before any worker search starts, its exclusive manifest exists continuously,
> and a read performed after its creation has established that the legacy
> driver's flushed START frontier is still earlier than that condition.

If legacy had passed the fresh-manifest guard before creation, its START flush
would already be visible in that second read, and the worker would abort. If
legacy reaches the condition later, the continuously present reservation makes
its fresh guard stop. These are the only two orderings under the stated source
and log semantics. A worker crash can leave pending work, but cannot leave the
legacy process suspended.

This proof requires the log to be the legacy process's **actual direct regular
stdout file**. A regular file populated by a delayed pipe, `tee`, PTY collector
or asynchronous forwarding service is insufficient. No independent process may
truncate the log or write fake START records. Keep worker logs separate. The
real backend must verify the binding; an assumed path or a plausible last line
does not earn it. Source hashes alone cannot prove the live file-descriptor
binding.

`SIGSTOP`/`SIGCONT` is unnecessary under this protocol. A `finally` block does not
recover from coordinator `SIGKILL` between STOP and CONT, and ordinary PID
check-then-signal has a reuse race without a stable process handle. Those risks
are avoided here by making no process-control calls.

## A reservation deliberately stops the original loop

When legacy reaches the first reserved prepared condition it will raise at its
fresh-manifest guard. That is an expected orchestration stop, not a method's
native failure. The original progress rows will then omit future conditions.
They **must not** become the completion denominator.

Before any reservation, prepare one preserving coordinator that owns the entire
unchanged manifest SHA256
`0dc9e3975d5451d3e21b2f474b7218ecf4d82e110e197fab95630bb13c22579b`,
modern configuration SHA256
`9d0d52456f2beab72cb007f5023663711aefa5eb55859ea845e4739190ceea2d`,
and VC configuration SHA256
`8b9b75de3fc0df376f867e2319d75b30b3f1f2b66a270f1a1af3e80c44454936`.
That coordinator must:

- Serialize reservation decisions for the two extra workers and keep at most
  two new workers active, with one native search per worker at CPU2.
- Reconcile all36 conditions from their raw modern/VC receipts, identities,
  source hashes and artifact pins. Fully terminal failures remain terminal.
  A status label alone does not establish a verified terminal case.
- Let active legacy and owned worker attempts finish. Never adopt or rerun a
  nonempty partial manifest, interrupted attempt or unknown live writer.
- After the expected legacy guard stop, continue scheduling every untouched
  condition, including deterministic preparation-failure records, through the
  same frozen execution contract. Any unexpected legacy exit or ownership
  ambiguity is a fail-stop requiring explicit inspection, not a retry.
- Preserve every unscheduled condition as pending and retain blocked conditions
  explicitly. Do not invent failure receipts for work that has not started.
- Claim complete only when all216 eligible slots are independently verified
  terminal. Native attempts, failures and unstarted slots remain distinct.

Thus no unstarted work disappears when the old loop stops. For example, five
verified early conditions plus two verified reverse-queue conditions account
for42 slots; the remaining29 conditions, 174 slots, are still pending. An owned
empty reservation is not a completed run. This recovery/coordinator work is a
required part of the optimization, not a follow-up that can be deferred until
after creating reservations.

## Synthetic evidence and recommendation

`aav_benchmark/paired_worksharing_model.py` enumerates all70 interleavings of the
four relevant legacy events and four worker events. The postclaim barrier
prevents a duplicate search in every model ordering. Removing that barrier
produces a concrete counterexample: legacy checks the missing manifest, then
the worker reserves and starts, then legacy writes its manifest and starts too.
The tests also exercise incomplete START prints, log rebinding/regression,
ownership loss, preparation-failure branches, future-only selection, worker
budget cleanup and preservation of the216-slot denominator.

Two controls invoke the unchanged modern engine only with mocked native
subprocesses in temporary directories: a fresh call rejects an existing empty
manifest without a search, and an owned-empty resumed call produces exactly
the three scheduled mock seed jobs. No actual native program runs.

```sh
python -m pytest -q tests/test_paired_worksharing_model.py
```

All12 tests pass. This is a proof of the stated atomic-event model and frozen
engine guard, not verification of a live stdout binding or a production
coordinator. The acceleration is worthwhile only if the remaining paired work
outlasts implementation/review of that coordinator. Keep the existing serial
driver running while this remains a proposal; enable nothing until the old156
gate and the complete preserving coordinator are both ready.
