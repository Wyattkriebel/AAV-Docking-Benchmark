# Field-specific source replay roundoff

The public twelve-reference replay allows an absolute **1e-12 Å** difference
only for the recomputed `native_bond_graph[*].distance_A` scalar. This field is
the Euclidean distance derived from pinned source atom coordinates. Deposited
`struct_conn_distance_A`, atom coordinates and identities, chemical connectivity,
source SHA pins, selected masks, assembly transforms and eligibility decisions
remain exact. The existing density tolerance and its exact threshold
classifications are unchanged.

The motivating observation is the [PR7 Linux CI failure](https://github.com/AlecKriebel/AAV-Docking-Benchmark/actions/runs/41689643039):
five source-replay tests stopped at the same derived bond distance,
1.4351466823986994 Å versus frozen 1.4351466823986996 Å, a difference of
2.22e-16 Å. Other tests included all 21 ring-initialization checks passing.
The previous source PR had passed Linux CI. Those observations do not establish
the cause of the later low-bit difference. No additional field receives a
tolerance merely because another difference might appear later.

The original atomic/reference helpers, ledgers, source coordinates and raw
archives remain unchanged. The public wrapper uses a scoped comparison hook;
it restores the original comparator after success or failure and does not relax
another thread's comparisons. A report records the configured tolerance, number
of compared and unequal fields, and maximum actual absolute bond-distance
difference separately from the density differences.

The default public profile is explicit in the API and command-line help. Request
strict equality with zero:

```sh
python scripts/verify_public_source_cohort.py --bond-distance-tolerance-A 0
```

For the public cross-platform profile:

```sh
python scripts/verify_public_source_cohort.py --bond-distance-tolerance-A 1e-12
```

Negative or nonfinite tolerances are rejected. Controls accept the observed
one-ULP derived difference, reject differences above the allowance, and reject
even one-ULP changes in raw coordinates, deposited distances or unrelated fields.
Graph-identity, decision, source-pin and density-classification controls remain
active. Private mutated-ledger fixtures use the same scoped source comparison
so that their intentional atom-loss and wrong-assembly checks still execute.

This is a source-validation portability change on the same frozen data. It does
not change any pose threshold, experimental qualification, native execution,
biological denominator or accuracy result. A downstream analysis plan that pins
the previous public-wrapper hash requires an explicit source-validation amendment
before binding a result to the revised wrapper.
