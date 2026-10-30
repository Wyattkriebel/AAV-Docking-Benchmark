# Implemented scope and remaining comparisons

The current product is an executable development diagnostic for **AutoDock Vina
1.2.7**, with structural and assay-reference evidence kept separate. It does not
complete the confirmatory design in `PROTOCOL.md`.

| Component | Current evidence | Remaining scientific requirement |
|---|---|---|
| AAV9 structure handling | Archived 7MTZ/7MT0, declared 60-mer, complete neighborhoods, independently solved apo, monomer controls, fixed native atom/contact maps | Independent local-map/density adjudication; local-interface alignment sensitivity; preparation/water sensitivity |
| Vina pose search | Three frozen seeds, two CPUs/run, oracle-site information, retained raw poses and full failure ledger | Blind/predicted-site search; an untouched, independently qualified cohort; convergence and time-budget comparisons |
| Linked glycans | Seven selected external examples, including one chemical-reference failure; six supported preparations | Local structural quality/assembly qualification; broader chemistry, including sulfation; independent family/fold splits |
| Native-pose scoring | Fixed experimental heavy atoms, explicitly modeled hydrogens, same receptor/box/score | Native controls are diagnostics, never docking success or proof about internal candidate sampling |
| Pocket methods | Historical fpocket/P2Rank files retained with unknown execution provenance | New versioned full-surface runs, symmetry-aware site evaluation and credible surface-patch baselines |
| Other docking/scoring methods | No new Vina-Carb, HADDOCK3, GNINA or GlycoTorch Vina execution | Versioned installation, released-example reproduction, supported chemistry, comparable information and compute budgets; learned-method training-overlap audit |
| Experimental affinity ranking | Three measured Gal3 binders from one ITC context; both reducing-end forms; size baseline | Assayed full-length construct differs from docking domain; three ligands cannot establish generalization; size alone explains their rank |
| Tested-negative specificity | Candidate public sources and their conflicts recorded in `NEXT_EVIDENCE.md` | Recover and qualify exact tested chemicals, presentation, context, uncertainty and labels before freezing a new panel |
| Uncertainty | Tested family-aware aggregation/bootstrap utilities | No family-level reliability interval is earned by these development inputs; source family names are not independently validated splits |

The six successfully prepared nonviral examples are not six established
independent families: 1PWB and 2VUZ share a broad C-type lectin fold. 2VUZ is a
resolved seven-residue fragment, not whole-nonasaccharide atomic truth. The
1QFO mismatch is retained in the seven-case eligible denominator. These facts
prevent interpreting successful execution as completion of the protocol's
independent-family and reference-quality gates.

All structural results use a dry, rigid receptor model. Retained calcium is a
fixed +2 Vina atom, without a specialized coordination or polarization model.
The neighborhood and box are reference-derived. No blind-site, unbound-input
generalization, glycan classifier, dissociation-constant prediction or universal
docking-reliability claim follows from this implementation milestone.

The remaining closure sequence is to independently qualify references and
assays, freeze family/chemical splits and estimands, reproduce each added
workflow, execute equally informed comparisons, retain all eligible failures,
and report held-out family-level uncertainty. Parameter changes made during
development must not be retroactively described as a frozen confirmatory run.
