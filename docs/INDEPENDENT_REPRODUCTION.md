# Independent AAV9 replay

On 28 October 2026 (America/Los_Angeles), a separate agent reran the frozen
`holo_complete` AAV9 diagnostic with Vina 1.2.7 and seed 104729. It invoked the
official executable directly and evaluated the returned PDBQT using a separate
parser and explicit RMSD/contact arithmetic. The reproduction does not import
the repository's engine or metrics modules.

The run reused the primary prepared receptor, independently embedded galactose,
reference and executable after checking their hashes against the primary run.
It changed only the output destination. The separate
[configuration](../results/independent_reproduction/aav9_one_seed/config.json)
contains seed 104729; the primary three-seed configuration and 12-run AAV9
denominator were preserved. No large inputs were copied.

| Endpoint | Independent result | Difference from primary |
|---|---:|---:|
| Returned poses | 20 | 0 |
| Top-1 Vina score (kcal/mol) | −4.865 | 0 |
| Top-1 fixed-frame ligand RMSD (Å) | 11.490781029445593 | 0 |
| Top-5 minimum RMSD (Å) | 8.473830646958517 | 0 |
| Best retained RMSD (Å) | 3.333033140649321 | 0 |
| Rank of best retained pose | 13 | 0 |
| Top-1 native-contact recall / precision | 0 / 0 | 0 / 0 |

All 20 ranks had identical scores, atom coordinates and independently computed
RMSDs. The reproduced pose file also has the same SHA256 digest as the primary
file: `2e2c5d7857f757f5ead40f44bd7b859a7711ae0bca8021c54e7dfbef6d956031`.
This exact replay verifies deterministic execution for this executable, host,
configuration and input set. It is not an additional biological replicate.

The independent controls established the following:

- The 12 native heavy-atom names and coordinates exactly match deposited
  7MTZ, author chain A, GAL 801. Duplicate or missing names fail evaluation.
- Native coordinates give 0 Å RMSD and all 35 native heavy-atom contact pairs
  at the strict `<4 Å` cutoff, with recall and precision both 1.
- Translating the native ligand by (3, 4, 0) Å gives 5 Å RMSD. No independent
  ligand superposition is used.
- Applying the same explicit proper rigid transform to the ligand, reference
  and receptor preserves RMSD and native contacts.
- Source-coordinate checksums match preparation. Applying the archived
  apo-to-holo rotation and translation to all 31,080 matched source Cα atoms
  reproduces the preparation RMSD of 0.6071591031339104 Å exactly. The rotation
  has determinant 1 and satisfies orthogonality within numerical precision.

The [machine-readable record](../results/independent_reproduction/aav9_one_seed/reproduction.json)
contains the command, hashes, timings, numerical differences and controls. The
[per-pose table](../results/independent_reproduction/aav9_one_seed/independent_pose_metrics.csv),
[raw output](../results/independent_reproduction/aav9_one_seed/holo_complete_104729.pdbqt)
and [execution log](../results/independent_reproduction/aav9_one_seed/holo_complete_104729.log)
are archived beside it. After preparing the frozen primary inputs, this check
is run with:

```sh
.venv/bin/python scripts/reproduce_aav9_one_seed.py
```

The script refuses to overwrite an existing replay artifact. Preserve or move
the replay directory before repeating it.

This is an independent implementation check on shared prepared inputs and the
same host. It does not independently adjudicate local density, reproduce all
preparation from scratch, execute six development systems, or test held-out
site recovery or assay specificity. It therefore supports the reproduction
component of closure work; it does not by itself pass the broader Gate 5 or
establish the repository's full scientific conclusion.
