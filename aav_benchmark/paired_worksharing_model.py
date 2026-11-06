"""Synthetic scheduling proof only: no signals, files, processes or native calls.

This model is not a live coordinator. Its predicates require independently
verified process/log ownership and source bindings. It cannot authorize or
perform a reservation, resume or native execution.
"""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass


PARENT_EVENTS = ("parent_START_flush", "parent_manifest_check", "parent_write_running", "parent_search")
WORKER_EVENTS = ("worker_O_EXCL", "worker_postclaim_START_read", "worker_owned_empty_check", "worker_search")


class SchedulingStop(ValueError):
    pass


@dataclass(frozen=True)
class LogFrontier:
    identity: tuple
    byte_length: int
    started_indices: tuple

    @property
    def index(self):
        return self.started_indices[-1] if self.started_indices else -1


def parse_frontier(raw, condition_ids, identity):
    lookup = {name: i for i, name in enumerate(condition_ids)}
    if len(lookup) != len(condition_ids):
        raise SchedulingStop("Duplicate condition identity")
    started = []
    for line in raw.splitlines(keepends=True):
        if line.endswith(b"\n") and line.startswith(b"START_MODERN "):
            name = line.decode("utf-8").rstrip("\r\n")[len("START_MODERN "):]
            if name not in lookup:
                raise SchedulingStop("Unplanned START condition")
            started.append(lookup[name])
    if started != list(range(len(started))):
        raise SchedulingStop("START records must be the one fixed forward-order prefix")
    return LogFrontier(tuple(identity), len(raw), tuple(started))


def postclaim_decision(before, after, index, *, prepared, exclusive_creation, owns_empty,
                       verified_regular_stdout_binding):
    """A pure decision, never a live claim or native-execution permission."""
    if not prepared or index <= before.index:
        raise SchedulingStop("Only a future prepared condition can use the legacy guard")
    if not verified_regular_stdout_binding:
        raise SchedulingStop("The exact regular stdout log must be independently bound")
    if not exclusive_creation:
        raise SchedulingStop("Existing output cannot be claimed or retried")
    if (after.identity != before.identity or after.byte_length < before.byte_length
            or after.started_indices[:len(before.started_indices)] != before.started_indices):
        raise SchedulingStop("Bound driver log changed identity or regressed")
    if after.index >= index:
        return "abort_and_retain_reservation"
    if not owns_empty:
        raise SchedulingStop("Owned empty manifest/receipt changed; no partial retries")
    return "barrier_holds_in_model"


def interleavings(left=PARENT_EVENTS, right=WORKER_EVENTS):
    if not left:
        yield tuple(right)
    elif not right:
        yield tuple(left)
    else:
        for rest in interleavings(left[1:], right):
            yield (left[0],) + rest
        for rest in interleavings(left, right[1:]):
            yield (right[0],) + rest


def simulate(schedule, *, postclaim_barrier=True):
    """Model one candidate and every ordering of legacy/worker atomic events."""
    if len(schedule) != 8 or set(schedule) != set(PARENT_EVENTS + WORKER_EVENTS):
        raise SchedulingStop("Only the eight declared atomic events are modeled")
    if tuple(e for e in schedule if e in PARENT_EVENTS) != PARENT_EVENTS:
        raise SchedulingStop("Parent event order changed")
    if tuple(e for e in schedule if e in WORKER_EVENTS) != WORKER_EVENTS:
        raise SchedulingStop("Worker event order changed")
    manifest = None
    parent_started = parent_blocked = worker_aborted = False
    claimed = parent_search = worker_search = False
    for event in schedule:
        if event == "parent_START_flush":
            parent_started = True
        elif event == "parent_manifest_check":
            parent_blocked = manifest is not None
        elif event == "parent_write_running" and not parent_blocked:
            manifest = "legacy_running"
        elif event == "parent_search" and not parent_blocked:
            parent_search = True
        elif event == "worker_O_EXCL":
            if manifest is None:
                manifest, claimed = "owned_empty", True
            else:
                worker_aborted = True
        elif event == "worker_postclaim_START_read" and claimed:
            if postclaim_barrier and parent_started:
                worker_aborted = True
        elif event == "worker_owned_empty_check" and claimed:
            if manifest != "owned_empty":
                worker_aborted = True
        elif event == "worker_search" and claimed and not worker_aborted:
            worker_search, manifest = True, "worker_running"
    return dict(parent_search=parent_search, worker_search=worker_search,
                duplicate_search=parent_search and worker_search,
                parent_guard_stopped=parent_blocked, reservation_created=claimed,
                any_process_paused=False)


def reconcile_schedule(condition_ids, verified_terminal_condition_ids, *, blocked_ids=()):
    """Preserve the full denominator even when the old sequential loop exits."""
    if len(condition_ids) != 36 or len(set(condition_ids)) != 36:
        raise SchedulingStop("The unchanged36-condition schedule is required")
    terminal = set(verified_terminal_condition_ids)
    blocked = set(blocked_ids)
    if (len(terminal) != len(verified_terminal_condition_ids) or terminal & blocked
            or not terminal <= set(condition_ids) or not blocked <= set(condition_ids)):
        raise SchedulingStop("Unknown or duplicate condition in recovery ledger")
    pending = [name for name in condition_ids if name not in terminal]
    return dict(eligible_method_seed_slots=216, verified_terminal_seed_slots=6 * len(terminal),
                pending_condition_ids=pending, blocked_condition_ids=sorted(blocked),
                complete=not pending and not blocked,
                requires_preserving_coordinator=bool(pending or blocked))


class SyntheticWorkerBudget:
    """Synthetic finally/CPU-budget check, deliberately without an OS backend."""
    def __init__(self, *, old156_verified_terminal, recovery_coordinator_ready):
        if old156_verified_terminal is not True or recovery_coordinator_ready is not True:
            raise SchedulingStop("Old156must be verified terminal and full recovery ready before any claim")
        self.active = 0

    @contextmanager
    def slot(self):
        if self.active >= 2:
            raise SchedulingStop("At most two additionalCPU2workers")
        self.active += 1
        try:
            yield
        finally:
            self.active -= 1
