# Implemented scope and remaining comparisons

The current product includes executed development diagnostics for **AutoDock Vina
1.2.7**, **fpocket 4.2.3** and **P2Rank 2.5.1**, with structural and assay-reference evidence kept separate. It does not
complete the confirmatory design in `PROTOCOL.md`.

| Component | Current evidence | Remaining scientific requirement |
|---|---|---|
| AAV9 structure handling | Archived 7MTZ/7MT0, declared 60-mer, complete neighborhoods, independently solved apo, monomer controls, fixed native atom/contact maps | Independent local-map/density adjudication; local-interface alignment sensitivity; preparation/water sensitivity |
| Vina pose search | Three frozen seeds, two CPUs/run, oracle-site information, retained raw poses and full failure ledger | Blind/predicted-site search; an untouched, independently qualified cohort; convergence and time-budget comparisons |
| Linked glycans | Seven selected external examples, including one chemical-reference failure; six supported preparations | Local structural quality/assembly qualification; broader chemistry, including sulfation; independent family/fold splits |
| Native-pose scoring | Fixed experimental heavy atoms, explicitly modeled hydrogens, same receptor/box/score | Native controls are diagnostics, never docking success or proof about internal candidate sampling |
| Pocket methods | New frozen whole-protein fpocket/P2Rank runs on nine conditions; source-only capsid symmetry, 100 area-weighted surface draws/condition; all runtime failures retained | Untouched families, broader assembly controls, additional pocket methods and comparable compute budgets |
| Other docking/scoring methods | Original Vina-Carb 1.0 tutorial reproduced with active carbohydrate terms and a validated fixed-pose scorer; HADDOCK3 2026.10.1 official guided test completed with audited topology, rebuilding and model counts | Paired new-cohort execution, supported chemistry, comparable information and compute budgets; GNINA/GlycoTorch remain unexecuted and require separate qualification |
| Experimental affinity ranking | Three measured Gal3 binders from one ITC context; both reducing-end forms; size baseline | Assayed full-length construct differs from docking domain; three ligands cannot establish generalization; size alone explains their rank |
| Tested-negative specificity | Six source-qualified Gal9C measured/NB epimer pairs; 36 frozen Vina runs; size/charge baselines; NB categorical at tested 2 mM | Untouched target families and chemical series, assay-construct sensitivity; threshold calibration and held-out sensitivity/specificity if a classifier is claimed |
| Uncertainty | Tested family-aware aggregation/bootstrap utilities | No family-level reliability interval is earned by these development inputs; source family names are not independently validated splits |

The six successfully prepared nonviral examples are not six established
independent families: 1PWB and 2VUZ share a broad C-type lectin fold. 2VUZ is a
resolved seven-residue fragment, not whole-nonasaccharide atomic truth. The
1QFO mismatch is retained in the seven-case eligible denominator. These facts
prevent interpreting successful execution as completion of the protocol's
independent-family and reference-quality gates.

All structural results use a dry, rigid receptor model. Retained calcium is a
fixed +2 Vina atom, without a specialized coordination or polarization model.
The original pose-diagnostic neighborhood and box are reference-derived. New
site prediction uses whole protein inputs, including independently solved AAV
apo; it does not yet supply blind-site docking pose recovery. The Gal9C box
uses a prior known site. No held-out-family generalization, calibrated glycan
classifier, dissociation-constant prediction or universal docking-reliability
claim follows from these development results.

The remaining closure sequence is to independently qualify references and
assays, freeze family/chemical splits and estimands, reproduce each added
workflow, execute equally informed comparisons, retain all eligible failures,
and report held-out family-level uncertainty. Parameter changes made during
development must not be retroactively described as a frozen confirmatory run.
