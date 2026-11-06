# Saved-product review provenance amendment

The initializer's original receipt pins 32 primary artifacts and a separate
derived independent-review file. A later standalone replay overwrote that
derived review after the receipt was created. Its author confirmed that the
replay changed provenance and description fields. The original review bytes
could not be recovered; an exact historical content comparison is unavailable.

The original receipt, its historical checksum and the current named review are
preserved. A [versioned copy](../results/ideal_ring_conformer_feasibility/independent_review_v2.json)
and an explicit [amendment](../results/ideal_ring_conformer_feasibility/review_pin_amendment.json)
record the replacement. The verifier reports
`historical_independent_review_pin_reproduced: false` separately from successful
verification of the versioned review. It does not treat the missing historical
bytes as reproduced or reconstruct a file to match their checksum.

All 32 primary initializer artifacts remain byte-identical, including the
algorithm, policy, seven selected conformers and preserved development failures.
The read-only verifier independently checks the complete chemical graphs,
freshly assigned three-dimensional stereochemistry, ring and clash gates,
candidate selection, mathematical controls and linkage-rotation controls.
It generates no conformer, docking pose or experimental accuracy result.
This later verification is a publication audit, not an external preregistration.

```sh
python scripts/verify_ring_initialization.py
```

The [original method description](INPUT_RING_GEOMETRY.md) and its own receipt
pin remain unchanged. This amendment corrects review provenance without
changing a frozen preparation or using docking performance to select inputs.
