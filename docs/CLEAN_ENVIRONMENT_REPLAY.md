# Clean-environment AAV9 replay

On 29 October 2026 (America/Los_Angeles), the complete AAV9 pipeline ran in a new
Python 3.11 virtual environment without inherited system packages. Installation
used `python -m pip install -e '.[test]'`. The isolated environment passed the
then-current 93-test suite; later controls and their final count are recorded
in the final implementation review and PR validation.

The replay restored all three AAV input files into a new raw-input directory
from the committed gzip archives, verifying their frozen SHA256 digests. It
prepared fresh receptors and ligand, performed all twelve configured docking
runs, and recomputed pose/contact metrics in a new output directory. It did
not restore or reuse the primary prepared receptor or ligand files.

Every returned pose file matches its primary-run counterpart byte for byte:
four apo/holo and neighborhood/monomer conditions, each with three seeds. The
[machine-readable comparison](../results/clean_environment/aav9/clean_environment_record.json)
contains all twelve primary/replay SHA256 pairs, software versions and runtime.
Every freshly generated preparation file also matches the primary archived
preparation byte for byte; all nineteen comparisons are retained in that record.
The primary input bundle therefore preserves the exact clean-replay input bytes.
Raw replay poses, logs, configurations and run manifests are archived alongside
that record. Wall time for source restoration through analysis was 169.156 s.

The environment used Python 3.11.8, NumPy 2.4.2, SciPy 1.17.0, Gemmi 0.7.5,
RDKit 2026.9.1 and Meeko 0.8.0. Pillow 13.0.0 and pytest 8.2.2 were installed
fresh. The executable was the same checksum-pinned official Vina 1.2.7 ARM64
binary on the same macOS host.

For another isolated replay, create a new environment, install the repository,
then run the documented `replay aav9` command with new `--prepared` and `--output`
directories. Input hashes are checked even when the original raw directory
already exists. `fetch` restores the archived source bytes when it does not.

This checks software installation, source restoration, preparation and seeded
execution on one host. It is not cross-platform reproduction, independent
structural-density adjudication, another biological observation, or a
confirmatory test cohort. It supplements the separate independent-parser
replay; neither check completes the repository's scientific purpose.
