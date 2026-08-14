# Validation Record

**Generated file — do not edit by hand.** Regenerate with:

```bash
python -m validation.scientific_method run
```

Every claim this project makes about itself is listed here with the
experiment that decides it and what that experiment actually measured.
A FALSIFIED row is a result, not a defect to hide: revise the claim with
`python -m validation.scientific_method revise <ID> --statement ...` and
run it again.

| Status | ID | Claim | Evidence |
|---|---|---|---|
| FALSIFIED | C11 | Despite non-convexity, annealing on E8 converges reliably -- GAS finds lower energy than uniform random sampling at equal evaluation budget. | GAS mean 0.280879 vs random 0.284600 (+1.3%), GAS wins 5/8 seeds |
| FALSIFIED | C15 | Empirical validation on known problems shows global recovery. | No validation study, dataset, benchmark or result table exists anywhere in the repository. The sentence cites work that was never done. |
| unsupported | C13 | With probability 1, GAS converges to a local minimum of E(x) on the E8 manifold. | Labelled a Theorem but supported only by a four-bullet sketch appealing to 'standard simulated annealing theory'. Those results require a cooling schedule of order 1/log(t); the implemented schedule is geometric, which does not satisfy their hypotheses. No proof is given here. |
| unsupported | C16 | Sphere packing, supply chain and energy system examples validate the framework against known optimal solutions. | No examples/ directory, no benchmark comparing against known optima. benchmarks/random_baseline.py compares GAS against random search only. |
| not testable here | C17 | Conventional equation architecture operates at 0.5 sigma quality, with a 94% defect rate across economics, AI training and corporate metrics. | No dataset, sample frame, population definition, coding protocol or inter-rater procedure is given for any of the percentages. The sigma levels are computed from figures that were not measured. |
| not testable here | C18 | Robots are 30x-100x more energy expensive than efficient human workers on a full lifecycle basis. | No citations, no lifecycle inventory, no system boundary and no source for any of the kWh/day figures. The ranges are asserted, not derived. |
| untested | C14 | For problems whose global optimum has rho_coset > rho_min, GAS finds the global minimum with high probability. | Explicitly labelled a conjecture, which is honest. Deciding it needs a problem family with known global optima; none exists in this repository yet. |
| PASS | C01 | The 240 generated vectors form the E8 root system: norm^2 = 2, integral inner products, no duplicates. | 240 roots, all norm^2=2.0, inner products in [-2.0, -1.0, 0.0, 1.0, 2.0], 0 duplicates |
| PASS | C02 | Each geometric energy term is a non-constant function of the state x. | worst term OctahedralEnergy: 300/300 distinct values over the sphere |
| PASS | C03 | No energy term degenerates to a constant over the sphere. | narrowest range: DodecahedralEnergy spans 0.003469 |
| PASS | C04 | Analytic gradients agree with central finite differences to better than 1e-5 relative error. | worst relative error 3.655e-09 (DodecahedralEnergy), tolerance 1e-5 |
| PASS | C05 | Every energy term has a nonzero gradient with respect to x. | smallest gradient norm across all terms and samples: 9.661e-03 (HexagonalEnergy) |
| PASS | C06 | The recorded energy changes over the course of an optimization run. | 273 distinct energies over 301 iterations |
| PASS | C07 | GAS reduces the geometric energy below its starting value. | mean improvement 14.1% over 8 seeds, worst 6.3% |
| PASS | C08 | The --all-terms flag changes the computation performed. | 3-term E=0.304758 (weights 3/3), 7-term E=0.535482 (weights 7/7) |
| PASS | C09 | The annealing schedule decays with the iteration index. | sigma_t 0.05000 -> 0.00204, T_t 0.004104 -> 0.000167 over 800 iterations |
| PASS | C10 | Runs are reproducible from a seed, independent of the global NumPy RNG state. | two runs under different global seeds: 0.275303544465 vs 0.275303544465 |
| PASS | C12 | Under the GAS update rule the expected energy decreases: E[E_{t+1}] <= E_t + O(sigma_t^2). | mean one-step dE = -0.003779 over 120 states, bound sigma^2 = 0.002500; 77% of steps were non-increasing |
| PASS | C19 | The meta-layer's geometric regularization term influences the decoded N-dimensional solution. | min \|\|grad R_geo\|\| = 3.032e-01; turning lambda_3 off moves the decoded solution by at least 7.166e-01 |

## Detail

### C11 — FALSIFIED

> Despite non-convexity, annealing on E8 converges reliably -- GAS finds lower energy than uniform random sampling at equal evaluation budget.

- **Asserted in:** README.md#why-this-matters
- **Kind:** empirical
- **Experiment:** `beats_random_search`
- **Revision:** 1
- **Measured:** GAS mean 0.280879 vs random 0.284600 (+1.3%), GAS wins 5/8 seeds
- **Note:** The headline performance claim. Requires GAS to win at least 6 of 8 seeds against the baseline.

### C15 — FALSIFIED

> Empirical validation on known problems shows global recovery.

- **Asserted in:** THEORY.md#9-convergence (9.3, listed as Evidence)
- **Kind:** provenance
- **Experiment:** `none available`
- **Revision:** 1
- **Measured:** No validation study, dataset, benchmark or result table exists anywhere in the repository. The sentence cites work that was never done.
- **Note:** Removed from THEORY.md. Recorded here so the removal is auditable rather than silent.

### C13 — unsupported

> With probability 1, GAS converges to a local minimum of E(x) on the E8 manifold.

- **Asserted in:** THEORY.md#9-convergence (9.2, stated as a Theorem)
- **Kind:** theoretical
- **Experiment:** `none available`
- **Revision:** 1
- **Note:** Labelled a Theorem but supported only by a four-bullet sketch appealing to 'standard simulated annealing theory'. Those results require a cooling schedule of order 1/log(t); the implemented schedule is geometric, which does not satisfy their hypotheses. No proof is given here.

### C16 — unsupported

> Sphere packing, supply chain and energy system examples validate the framework against known optimal solutions.

- **Asserted in:** README.md#examples
- **Kind:** provenance
- **Experiment:** `none available`
- **Revision:** 1
- **Measured:** No examples/ directory, no benchmark comparing against known optima. benchmarks/random_baseline.py compares GAS against random search only.
- **Note:** Rewritten in README.md as planned work rather than completed validation.

### C17 — not testable here

> Conventional equation architecture operates at 0.5 sigma quality, with a 94% defect rate across economics, AI training and corporate metrics.

- **Asserted in:** legacy/Six-Sigma.md
- **Kind:** empirical
- **Experiment:** `none available`
- **Revision:** 1
- **Measured:** No dataset, sample frame, population definition, coding protocol or inter-rater procedure is given for any of the percentages. The sigma levels are computed from figures that were not measured.
- **Note:** Moved to legacy/. To make this testable, define the population of 'equations' being sampled, the operational definition of a 'defect', and the coding procedure -- then revise this claim and register an experiment.

### C18 — not testable here

> Robots are 30x-100x more energy expensive than efficient human workers on a full lifecycle basis.

- **Asserted in:** legacy/papers/energy.md
- **Kind:** empirical
- **Experiment:** `none available`
- **Revision:** 1
- **Measured:** No citations, no lifecycle inventory, no system boundary and no source for any of the kWh/day figures. The ranges are asserted, not derived.
- **Note:** Moved to legacy/. To make this testable, cite a published lifecycle assessment, state the system boundary, and register an experiment recomputing the ratio from the cited inventory.

### C14 — untested

> For problems whose global optimum has rho_coset > rho_min, GAS finds the global minimum with high probability.

- **Asserted in:** THEORY.md#9-convergence (9.3, stated as a Conjecture)
- **Kind:** conjecture
- **Experiment:** `none available`
- **Revision:** 1
- **Note:** Explicitly labelled a conjecture, which is honest. Deciding it needs a problem family with known global optima; none exists in this repository yet.

### C01 — PASS

> The 240 generated vectors form the E8 root system: norm^2 = 2, integral inner products, no duplicates.

- **Asserted in:** THEORY.md#2-e8-lattice, gas/lattice.py
- **Kind:** empirical
- **Experiment:** `e8_root_system_valid`
- **Revision:** 1
- **Measured:** 240 roots, all norm^2=2.0, inner products in [-2.0, -1.0, 0.0, 1.0, 2.0], 0 duplicates

### C02 — PASS

> Each geometric energy term is a non-constant function of the state x.

- **Asserted in:** THEORY.md#4-seed-equations, gas/energy_terms.py
- **Kind:** empirical
- **Experiment:** `energy_varies_with_x`
- **Revision:** 2
- **Measured:** worst term OctahedralEnergy: 300/300 distinct values over the sphere
- **Note:** Revision 1 of the code failed this: terms read only the neighbour set, never x.
- **Superseded wordings:**
  - r1 (FALSIFIED): "Each geometric energy term is a non-constant function of the state x." — Claim retained unchanged; the implementation was corrected instead (gas/energy_terms.py revision 2).

### C03 — PASS

> No energy term degenerates to a constant over the sphere.

- **Asserted in:** gas/energy_terms.py
- **Kind:** empirical
- **Experiment:** `no_term_is_constant`
- **Revision:** 2
- **Measured:** narrowest range: DodecahedralEnergy spans 0.003469
- **Note:** Targets the specific revision-1 collapses: TetrahedralEnergy == 1/3 and GoldenEnergy == |1-phi|.
- **Superseded wordings:**
  - r1 (FALSIFIED): "No energy term degenerates to a constant over the sphere." — Implementation corrected: tetrahedral now scores x-to-neighbour cosines, golden now scores shell-distance ratios from x.

### C04 — PASS

> Analytic gradients agree with central finite differences to better than 1e-5 relative error.

- **Asserted in:** gas/energy_terms.py
- **Kind:** empirical
- **Experiment:** `gradients_match_finite_difference`
- **Revision:** 1
- **Measured:** worst relative error 3.655e-09 (DodecahedralEnergy), tolerance 1e-5
- **Note:** Standard code verification. Revision 1 had no analytic gradient to check.

### C05 — PASS

> Every energy term has a nonzero gradient with respect to x.

- **Asserted in:** gas/energy_terms.py
- **Kind:** empirical
- **Experiment:** `gradients_nonzero`
- **Revision:** 2
- **Measured:** smallest gradient norm across all terms and samples: 9.661e-03 (HexagonalEnergy)
- **Superseded wordings:**
  - r1 (FALSIFIED): "Every energy term has a nonzero gradient with respect to x." — Implementation corrected; gradients are now analytic.

### C06 — PASS

> The recorded energy changes over the course of an optimization run.

- **Asserted in:** gas/solver.py
- **Kind:** empirical
- **Experiment:** `solver_energy_changes`
- **Revision:** 2
- **Measured:** 273 distinct energies over 301 iterations
- **Superseded wordings:**
  - r1 (FALSIFIED): "The recorded energy changes over the course of an optimization run." — step() now re-evaluates the proposal under its own neighbourhood.

### C07 — PASS

> GAS reduces the geometric energy below its starting value.

- **Asserted in:** README.md#why-this-matters
- **Kind:** empirical
- **Experiment:** `solver_descends`
- **Revision:** 1
- **Measured:** mean improvement 14.1% over 8 seeds, worst 6.3%

### C08 — PASS

> The --all-terms flag changes the computation performed.

- **Asserted in:** gas/cli.py, README.md
- **Kind:** empirical
- **Experiment:** `all_terms_flag_has_effect`
- **Revision:** 2
- **Measured:** 3-term E=0.304758 (weights 3/3), 7-term E=0.535482 (weights 7/7)
- **Superseded wordings:**
  - r1 (FALSIFIED): "The --all-terms flag changes the computation performed." — _compute_weights now returns one weight per term and raises on a count mismatch.

### C09 — PASS

> The annealing schedule decays with the iteration index.

- **Asserted in:** THEORY.md#9-convergence, gas/solver.py
- **Kind:** empirical
- **Experiment:** `annealing_schedule_decays`
- **Revision:** 2
- **Measured:** sigma_t 0.05000 -> 0.00204, T_t 0.004104 -> 0.000167 over 800 iterations
- **Note:** THEORY.md 9.2 rests on this: its proof sketch asserts sigma_t, T_t -> 0.
- **Superseded wordings:**
  - r1 (FALSIFIED): "The annealing schedule decays with the iteration index." — A cooling factor exp(-t/tau_anneal) was added to eta_t, alpha_t, sigma_t and T_t.

### C10 — PASS

> Runs are reproducible from a seed, independent of the global NumPy RNG state.

- **Asserted in:** gas/solver.py, gas/cli.py
- **Kind:** empirical
- **Experiment:** `reproducible_under_seed`
- **Revision:** 2
- **Measured:** two runs under different global seeds: 0.275303544465 vs 0.275303544465
- **Superseded wordings:**
  - r1 (FALSIFIED): "Runs are reproducible from a seed." — The solver now takes an injectable numpy Generator.

### C12 — PASS

> Under the GAS update rule the expected energy decreases: E[E_{t+1}] <= E_t + O(sigma_t^2).

- **Asserted in:** THEORY.md#9-convergence (9.1, stated as a Lemma)
- **Kind:** theoretical
- **Experiment:** `expected_energy_decreases`
- **Revision:** 2
- **Measured:** mean one-step dE = -0.003779 over 120 states, bound sigma^2 = 0.002500; 77% of steps were non-increasing
- **Note:** Revision-1 premise (deterministic gradient descent) was false because gradients were identically zero; with revision-2 energy terms the premise holds, so the lemma is now empirically decidable. Registered experiment expected_energy_decreases.
- **Superseded wordings:**
  - r1 (UNTESTED): "Under the GAS update rule the expected energy decreases: E[E_{t+1}] <= E_t + O(sigma_t^2)." — Revision-1 premise (deterministic gradient descent) was false because gradients were identically zero; with revision-2 energy terms the premise holds, so the lemma is now empirically decidable. Registered experiment expected_energy_decreases.

### C19 — PASS

> The meta-layer's geometric regularization term influences the decoded N-dimensional solution.

- **Asserted in:** THEORY.md#8-meta-layer, meta_layer/decoder.py
- **Kind:** empirical
- **Experiment:** `decoder_geometric_regularizer_active`
- **Revision:** 1
- **Measured:** min ||grad R_geo|| = 3.032e-01; turning lambda_3 off moves the decoded solution by at least 7.166e-01
- **Note:** With revision-1 energy terms the R_geo gradient measured ~3e-11, so lambda_3 was inert and the 'geometric coherence' regularizer did nothing.
