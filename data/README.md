The two blank AAV mmCIF placeholders were removed. `python -m aav_benchmark fetch`
restores verified coordinate files under `data/raw/` from the committed gzip
snapshots in `data/archive/`, or downloads the frozen source if no snapshot is
available. SHA256 checksums refer to the original uncompressed bytes.

Legacy nonviral mmCIF files and trimmed PDBs are historical inputs, not newly
qualified atomic references. The active development sources and their exact
graphs are listed in `benchmarks/external_cohort.json` and fetched from the
published HADDOCK dataset at an immutable Git commit. Its license is preserved
alongside archived source inputs. Reuse does not reproduce HADDOCK performance.

Generated preparation files live in `data/processed/`; final archived preparation
bundles are retained with each committed result. RCSB/wwPDB data and external
dataset provenance retain their own attribution; the repository code license
does not relicense third-party scientific data or Vina.
