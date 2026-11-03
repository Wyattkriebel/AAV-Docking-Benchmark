# Paired guidance and preparation contrasts

This bounded postprocessor aggregates recovery indicators from the four globally released structural JSON reports. It adds main-engine oracle-minus-fpocket/P2Rank contrasts within each engine/profile, and repair-minus-strict preparation sensitivities within each engine/site or the separate HADDOCK protocol. It does not recalculate a pose metric, fit, contact, native score or experimental correspondence.

Protocol labels remain H1 true-site guidance, H2 candidate ranking, H3 assay specificity and H4 assembly. Repair is a preparation sensitivity, never H3. Existing HADDOCK oracle-minus-predicted intervals, H2 common-union coefficient intervals and H4 assembly intervals are inventoried and referenced without recomputation. HADDOCK's existing `conditional_family_statistics` means a purposive-cohort interpretation; its failures already remain in the eligible denominator. Its AIR patches are not equivalent to hard Vina boxes, and no cross-backend comparison is added.

The prospective input packet is fixed to:

```text
results/structural_pose_v3/global_reports/main_paired.json
results/structural_pose_v3/global_reports/haddock.json
results/structural_pose_v3/global_reports/assembly_control.json
results/structural_pose_v3/global_reports/fixed_union.json
results/structural_pose_v3/global_report_release/generation_receipt.json
```

Root supplies that output directory to the unchanged global guard and saves its successful generation JSON at the receipt path. The postprocessor verifies `all_four_reports_generated`, all four ready readers, their fixed terminal denominators, the frozen global contract and authorization, exact output directory, report roles, complete report headers and metric-code inventories. The source numeric inputs consist only of these four reports; plans, freezes and code pins are metadata dependencies. The old156 development experiment remains with its existing reader.

The frozen global generation receipt does not contain report-file hashes. A separate helper therefore creates `results/structural_pose_v3/global_report_release/released_reports_binding.json` **after generation**. It binds the exact four filenames, roles, SHA-256 digests and byte counts, plus the generation-receipt digest and frozen global metadata. These hashes are explicitly post-generation bindings, not contemporaneous producer checksums or independent authentication of the earlier production event. Aggregation requires the binding's explicit SHA-256 and rechecks every bound file. The four-report packet remains unchanged.

All six cases, both profiles, three site arms and three seeds must be present exactly once. Main engines and HADDOCK remain separate. Five full-ligand cases and the single 11-atom modeled-core case stay separate. Unknown, missing or duplicated groups/seeds/engines fail validation, as do non-Boolean indicators, inconsistent rank budgets or failed outcomes with success. Partial HADDOCK physical models do not turn failed execution into a completed scored pool. H2's identical union and best-retained invariants are checked while its statistics are referenced.

There are 288 new main H1 strata and 324 preparation-sensitivity strata (216 main, 108 HADDOCK), covering global/local frames, Top1/Top5/best-retained and 1/2/3 Å. Each eligible paired estimate uses the unchanged whole-source-family bootstrap: seeds within case, cases within family, 10,000 requested replicates, seed20261103 and 95% confidence. Failed/unsupported outcomes remain zero. The single core family has no between-family interval. The intervals describe this purposive sample, not population or model-training-unseen reliability.

Completion summaries are separate: eligible pairs, completed A/B, both completed, only A, only B and neither; arm-specific completed-only rates and per-case completion counts. Completed-only rates average completed seeds within cases and cases within source families, so their contributing cohorts and denominators may differ between arms. No conditional difference or confidence interval is invented.

Before any structural accuracy is generated, root freezes `benchmarks/structural_contrast_evaluation.json` with `frozen_before_any_structural_accuracy=true` and an independent receipt containing its `aggregation_plan_sha256`, `frozen_before_any_structural_accuracy=true` and `new_structural_accuracy_consulted=false`. Then, only after the successful global generation:

```sh
python scripts/regenerate_structural_contrast_report.py bind --root . --plan FINAL_PLAN --freeze FREEZE
python scripts/regenerate_structural_contrast_report.py aggregate --root . --plan FINAL_PLAN --freeze FREEZE --binding-sha256 BINDING_SHA --output results/structural_pose_v3/structural_contrasts.json
```

The binding command prints the binding descriptor. The aggregation command writes a fresh result or verifies an identical prior artifact. It refuses to replace different existing output or write inside the four-report packet. Implementation and tests use source inspection, synthetic indicators and synthetic receipt packets only; no actual released metrics are consumed while preparing this proposal. The CLI and aggregation require only the Python standard library.
