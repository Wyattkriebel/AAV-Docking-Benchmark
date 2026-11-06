# Independent seed execution

Long cases can execute independent seeds concurrently while each run keeps its
frozen two-CPU search budget. This changes scheduling, not ligand inputs,
information, search settings, seeds or the scientific denominator. Wall-time
sums are not elapsed study duration when runs overlap.

After preparing a case, separate worker processes can run individual seeds:

```sh
python -m aav_benchmark dock 2VUZ --only-seed 104729 --output results/workers/2VUZ/104729
python -m aav_benchmark dock 2VUZ --only-seed 130363 --output results/workers/2VUZ/130363
python -m aav_benchmark dock 2VUZ --only-seed 155921 --output results/workers/2VUZ/155921
```

Each worker copies the full frozen three-seed configuration and records its
actual command, input hashes, poses, logs and resource use. `analyze` refuses a
worker's incomplete receptor-by-seed denominator. Once every worker has stopped:

```sh
python -m aav_benchmark merge 2VUZ --output results/replay/2VUZ --shard results/workers/2VUZ/104729 --shard results/workers/2VUZ/130363 --shard results/workers/2VUZ/155921
python -m aav_benchmark analyze 2VUZ --output results/replay/2VUZ
```

The merger requires exact complete coverage, one terminal record per eligible
identity, matching prepared inputs/configuration/executable fingerprints,
complete atom maps and matching pose hashes. It archives the source manifests,
preserves original commands and timings, and backs up changed existing artifacts.
Failed runs retain their places in the denominator. Duplicate or partial
coverage, changed inputs and live worker records fail visibly.

For a scheduler handoff, `--include-stopped-baseline` includes an existing output
manifest after its writer has been stopped. Interrupted baseline entries become
explicit prior attempts when a separately executed seed replaces them. Unknown
interruption timings stay unknown. They are never filled with zero or counted
as successful executions. The actual 2VUZ handoff is recorded separately with
the preserved completed first seed and any cancelled duplicate work.
