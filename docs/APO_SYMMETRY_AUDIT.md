# Source-only AAV apo symmetry audit

Audit and amendment date: 30 October 2026. Approximate site-class deduplication is
supported for the existing expanded AAV9 apo model, 7MT0, under an explicitly
amended **0.1 Å maximum C-alpha correspondence guard**, applied uniformly to
7MT0 and holo 7MTZ. The scientific 4 Å center-to-native endpoint and 4 Å candidate
merging distance remain unchanged, as do the 3/5 Å endpoint sensitivities.

**NoOutcomesUsed:** the recommendation and its adoption used only the pinned
receptor source coordinates and residue identities. The auditor did not inspect
predictor centers, candidate ranks, predictor recovery, or random-baseline
recovery. The coordinating agent accepted the guard amendment before any site
endpoint scoring and records the authoritative amendment and run receipt
separately. The original configuration and earlier failed preparation/baseline
records are retained. The apo observation remains in the experiment denominator.

Full source checksums, every operator's metrics and transforms, chain maps,
residual quantiles, and the separate later diagnostic appear in
`benchmarks/sources/apo_symmetry_audit.json`. Its
`source_only_decision_payload_sha256` hashes the source-only decision, method,
and source audits, excluding the later native-glycan diagnostic.

## Why the original guard failed

The initial procedure fits one receptor subunit onto each of the other subunits
using matching C-alpha residue identities. That pairwise fit proposes the
correspondence of every C-alpha in the entire capsid. An identity-preserving,
bijective whole-capsid fit then refines the operator; the refinement changes no
source coordinates or candidate centers.

All 60 subunits within each source have the same 518 named C-alpha identities.
The largest single-subunit fit RMSD is approximately 0.00051 Å for apo and
0.00053 Å for holo, close to the precision of the deposited coordinates. The
subunits within a given deposited state therefore show essentially rigid copies
of the same local fold. Fitting just one subunit can nevertheless amplify a small
angular discrepancy across the full capsid diameter. The global refinement
reduces that amplification, while residual approximate consistency of the
expanded subunit placements remains.

| Source-only quantity | Apo 7MT0 | Holo 7MTZ |
|---|---:|---:|
| Protein C-alpha atoms | 31,080 | 31,080 |
| Canonical protein heavy atoms | 247,860 | 247,860 |
| Maximum pairwise subunit C-alpha RMSD, Å | 0.0005104 | 0.0005262 |
| Maximum globally refined whole-capsid C-alpha displacement, Å | 0.0733156 | 0.0299343 |
| Maximum whole-capsid C-alpha RMSD across operators, Å | 0.0091231 | 0.0062370 |
| Operators exceeding original 0.05 Å maximum guard | 59/60 | 0/60 |
| Operators satisfying amended 0.1 Å maximum guard | 60/60 | 60/60 |
| Maximum canonical heavy-atom displacement under the same operators, Å | 0.0836416 | 0.0305957 |
| Maximum whole-capsid heavy-atom RMSD across operators, Å | 0.0092227 | 0.0062191 |

For the worst apo operator, fewer than 0.58% of C-alpha atoms exceed the old
0.05 Å guard. The full per-operator distributions and affected atom identities
are preserved rather than omitting difficult atoms or fitting a selected core.
This pattern supports approximate NCS placement/rotation consistency as the
issue. It does not support attributing the failure to meaningful local apo-fold
differences between individual symmetry copies. Real **apo-versus-holo**
conformational differences are a separate comparison and are not used to derive
the guard.

## Whole-capsid correspondence and the accepted scope

For both source models, all 60 refined operators preserve the residue identity of
every C-alpha and give bijections over the whole receptor. Each operator also
maps each entire chain to one entire chain. The 60 resulting chain permutations
are distinct, include identity, and close under composition. Their orders are
one identity, 15 operations of order two, 20 of order three, and 24 of order five:
the expected icosahedral permutation pattern. This is stronger than assuming that
60 similar-looking or homologous chains establish a capsid symmetry.

The additional heavy-atom diagnostic uses these same fixed operators and named
chain/atom correspondences. It does not refit heavy atoms. Every canonical
protein heavy atom also lies within 0.1 Å of its mapped source counterpart in
both states. No ligand, native site, mutagenesis residue, pocket prediction, or
algorithm outcome participates in deriving these operators or checks.

The 0.1 Å guard is an operational allowance for approximate receptor geometry,
not an accuracy threshold for a docking or pocket method. It is 2.5% of the
already frozen 4 Å merging distance, and one tenth of the 1 Å spacing of the
prespecified endpoint sensitivities. It allows source correspondence noise at a
small scale relative to this site-localization experiment while retaining a
strict whole-receptor maximum, residue identities, and bijectivity. It is not a
confidence interval or proof of exact atomic symmetry. Classification extremely
close to a merging boundary remains an approximate geometric decision.

The same bound applies to both holo and apo; the holo operators do not
numerically change. Scientific distance thresholds, candidate ranks, receptor
coordinates, method footprints, and native reference eligibility do not change.
The source-only implementation's single tolerance parameter may also guard the
initial subunit RMSD; its actual fits are far below either 0.05 or 0.1 Å. The
runtime receipt must state that parameter behavior accurately.

If a receptor still fails this amended, identity-preserving whole-capsid guard,
its symmetry-collapsed Top-5 endpoint must be marked unavailable or failed under
the frozen reporting rules. Such a case must not be silently dropped from the
experiment denominator. The successful source audit here permits approximate
site-class deduplication for 7MT0; it does not independently qualify an atomic
ligand pose.

## Later native-orbit diagnostic

Only after the source-only recommendation was accepted, a separate diagnostic
compared the fixed apo operators with the 60 deposited holo GAL copies. First,
a proper apo-to-holo fit used all 31,080 common protein C-alpha identities; its
RMSD was 0.607159 Å. This is the receptor-frame alignment, not a ligand fit.

The representative holo GAL was transformed into the apo frame, moved by each
fixed apo operator, and transformed back into the holo frame. Corresponding atom
names were compared against every deposited holo GAL copy. The 60 images match
60 distinct native copies. The largest atom displacement was 0.0197393 Å and the
largest corresponding 12-atom RMSD was 0.0192741 Å.

Those numbers are a post-adoption source-reproduction check. They did not select
operators, set the 0.1 Å guard, choose a predicted pocket, establish atomic pose
quality, or measure predictor accuracy. They are explicitly separated from the
source-only acceptance data in the machine-readable audit.
