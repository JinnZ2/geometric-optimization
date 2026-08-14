# Legacy

Superseded material, kept so that the record of what changed and why stays
auditable. Nothing here is imported by the package or exercised by CI.

Each item below is cross-referenced to the claim it relates to in
[`../validation/claims.json`](../validation/claims.json); current verdicts are
in [`../VALIDATION.md`](../VALIDATION.md).

## `gas/energy_terms.py`, `gas/solver.py` — revision 1

The original Geometric Annealing Solver. Superseded because an audit found the
optimization loop did not function:

| Defect | Measured | Claim |
|---|---|---|
| Energy terms never read the state `x` | `E(x)` identical for `x` and an unrelated vector, across all 7 terms | C02 |
| `TetrahedralEnergy` was the constant `1/3` | λ_min of a rank-≤8 Gram matrix of 24 vectors is always 0 | C03 |
| `GoldenEnergy` was the constant `\|1-φ\|` | it ratioed root *norms*, and every E8 root has norm √2 | C03 |
| Gradients were identically zero | finite differences perturbed `x` while holding the neighbour set fixed | C05 |
| Proposals scored against the *current* point's neighbours | `ΔE ≡ 0`; 1000 iterations gave 1 distinct energy | C06 |
| `--all-terms` was a no-op | `_compute_weights` returned 3 weights regardless of term count, so `zip()` dropped terms 4–7 | C08 |
| The annealing schedule never annealed | `σ_t` and `T_t` depended only on ρ, never on `t` | C09 |
| Sampling used the global `np.random` state | seeding was process-global | C10 |

The tests in `tests/test_energy_terms.py::TestFalsification` and
`tests/test_solver.py::TestFalsification` are written to fail against this code.
That is deliberate: a regression test that cannot fail is not evidence.

## `Six-Sigma.md`

Applies Six Sigma vocabulary to "equation architecture", reporting figures such
as a 94% defect rate and a 0.5σ baseline quality level.

Moved here because no dataset, population definition, sampling frame, defect
coding protocol or inter-rater procedure is given for any of the percentages —
the sigma levels are computed from numbers that were not measured. The document
also contains a section framing rejection of its argument as an admission of bad
faith, which is rhetoric rather than quality analysis.

Tracked as claim **C17** (`UNFALSIFIABLE_HERE`). To bring it back: define the
population of equations being sampled, the operational definition of a defect,
and the coding procedure, then revise the claim and register an experiment.

## `papers/energy.md`

A thermodynamic argument that efficient human workers outperform robots by
30×–100× on a full lifecycle basis.

Moved out of `papers/` because that directory implies a peer-reviewed standing
the document does not have. No citations, no system boundary, and no lifecycle
inventory are given for any of the kWh/day figures; the ranges are asserted
rather than derived.

Tracked as claim **C18** (`UNFALSIFIABLE_HERE`). To bring it back: cite a
published lifecycle assessment, state the system boundary, and register an
experiment that recomputes the ratio from the cited inventory.

---

Neither essay is claimed here to be *wrong*. Both are unfalsifiable as written,
which is a different and more fixable problem: there is no stated method by
which anyone — including their authors — could check them.
