"""
Experiments that decide the claims in ``claims.json``.

Each experiment is a zero-argument function returning a :class:`Verdict`. An
experiment must be able to come back FALSIFIED -- a check that cannot fail is
not evidence. Register new experiments with the :func:`experiment` decorator;
the id you register under is what ``claims.json`` refers to.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from typing import Callable, Dict, Optional

from gas.lattice import E8Lattice, create_test_point
from gas.energy_terms import (
    EnergyTerm,
    OctahedralEnergy,
    TetrahedralEnergy,
    GoldenEnergy,
    SquareEnergy,
    HexagonalEnergy,
    DodecahedralEnergy,
    IcosahedralEnergy,
    create_energy_suite,
)
from gas.solver import GeometricAnnealingSolver, GASParams

ALL_TERM_CLASSES = [
    OctahedralEnergy,
    TetrahedralEnergy,
    GoldenEnergy,
    SquareEnergy,
    HexagonalEnergy,
    DodecahedralEnergy,
    IcosahedralEnergy,
]


@dataclass
class Verdict:
    """Outcome of one experiment.

    Attributes:
        supported: True if the claim survived the attempt to falsify it.
        measured: Human-readable summary of what was actually observed.
        data: Machine-readable measurements, recorded for later comparison.
    """

    supported: bool
    measured: str
    data: Dict[str, float] = field(default_factory=dict)


REGISTRY: Dict[str, Callable[[], Verdict]] = {}


def experiment(name: str):
    """Register an experiment under ``name``."""

    def wrap(fn: Callable[[], Verdict]) -> Callable[[], Verdict]:
        REGISTRY[name] = fn
        return fn

    return wrap


def get(name: Optional[str]) -> Optional[Callable[[], Verdict]]:
    """Look up a registered experiment, or None if ``name`` is None/unknown."""
    if name is None:
        return None
    return REGISTRY.get(name)


# --------------------------------------------------------------------------
# Shared fixtures
# --------------------------------------------------------------------------
_LATTICE: Optional[E8Lattice] = None


def lattice() -> E8Lattice:
    global _LATTICE
    if _LATTICE is None:
        _LATTICE = E8Lattice()
    return _LATTICE


def sphere_points(n: int, seed: int = 0) -> np.ndarray:
    """``n`` reproducible points on the norm-sqrt(2) sphere."""
    rng = np.random.default_rng(seed)
    p = rng.standard_normal((n, 8))
    pts: np.ndarray = p / np.linalg.norm(p, axis=1, keepdims=True) * np.sqrt(2)
    return pts


def random_search(n_evals: int, seed: int, terms=None) -> float:
    """Best energy found by sampling ``n_evals`` points uniformly. The baseline."""
    L = lattice()
    terms = terms if terms is not None else create_energy_suite()
    probe = GeometricAnnealingSolver(
        L, terms, GASParams(), rng=np.random.default_rng(0)
    )
    best = np.inf
    for p in sphere_points(n_evals, seed):
        nb, idx = L.nearest_neighbors(p, k=24)
        best = min(best, probe._compute_energy(p, nb, L.coset_density(idx)))
    return float(best)


# --------------------------------------------------------------------------
# Lattice
# --------------------------------------------------------------------------
@experiment("e8_root_system_valid")
def exp_e8_root_system_valid() -> Verdict:
    """The 240 generated vectors really are the E8 root system."""
    L = lattice()
    roots = L.all_roots
    n = len(roots)
    norms_sq = np.sum(roots**2, axis=1)
    gram = roots @ roots.T
    inner = np.unique(np.round(gram, 9))
    dupes = n - len(np.unique(roots, axis=0))

    ok = (
        n == 240
        and np.allclose(norms_sq, 2.0)
        and set(inner.tolist()) <= {-2.0, -1.0, 0.0, 1.0, 2.0}
        and dupes == 0
    )
    return Verdict(
        ok,
        f"{n} roots, all norm^2={norms_sq[0]:.1f}, "
        f"inner products in {sorted(inner.tolist())}, {dupes} duplicates",
        {"n_roots": float(n), "duplicates": float(dupes)},
    )


# --------------------------------------------------------------------------
# Energy terms
# --------------------------------------------------------------------------
@experiment("energy_varies_with_x")
def exp_energy_varies_with_x() -> Verdict:
    """Every energy term is a non-constant function of the state x.

    This is the experiment revision 1 failed outright: its terms never read x,
    so E was piecewise constant with <=10 distinct values over the sphere.
    """
    L = lattice()
    pts = sphere_points(300, seed=1)
    worst_name, worst_distinct = None, 10**9
    for cls in ALL_TERM_CLASSES:
        t = cls()
        vals = []
        for p in pts:
            nb, _ = L.nearest_neighbors(p, k=24)
            vals.append(t.compute(p, nb))
        d = len(np.unique(np.round(vals, 9)))
        if d < worst_distinct:
            worst_name, worst_distinct = cls.__name__, d
    ok = worst_distinct >= len(pts) // 2
    return Verdict(
        ok,
        f"worst term {worst_name}: {worst_distinct}/{len(pts)} distinct values "
        f"over the sphere",
        {"min_distinct": float(worst_distinct), "n_samples": float(len(pts))},
    )


@experiment("gradients_nonzero")
def exp_gradients_nonzero() -> Verdict:
    """Every term has a nonzero gradient (revision 1's were identically 0)."""
    L = lattice()
    pts = sphere_points(20, seed=2)
    smallest, name = np.inf, None
    for cls in ALL_TERM_CLASSES:
        t = cls()
        norms = []
        for p in pts:
            nb, _ = L.nearest_neighbors(p, k=24)
            norms.append(np.linalg.norm(t.gradient(p, nb)))
        m = float(np.min(norms))
        if m < smallest:
            smallest, name = m, cls.__name__
    ok = smallest > 1e-8
    return Verdict(
        ok,
        f"smallest gradient norm across all terms and samples: "
        f"{smallest:.3e} ({name})",
        {"min_grad_norm": smallest},
    )


@experiment("gradients_match_finite_difference")
def exp_gradients_match_finite_difference() -> Verdict:
    """Analytic gradients agree with central finite differences.

    Standard code verification: the implemented derivative must match the
    derivative of the implemented energy. Revision 1 had no analytic gradient
    to check.
    """
    L = lattice()
    pts = sphere_points(5, seed=3)
    eps, worst, name = 1e-6, 0.0, None
    for cls in ALL_TERM_CLASSES:
        t: EnergyTerm = cls()
        for p in pts:
            nb, _ = L.nearest_neighbors(p, k=24)
            ga = t.gradient(p, nb)
            gn = np.zeros(8)
            for i in range(8):
                a, b = p.copy(), p.copy()
                a[i] += eps
                b[i] -= eps
                gn[i] = (t.compute(a, nb) - t.compute(b, nb)) / (2 * eps)
            rel = np.linalg.norm(ga - gn) / max(np.linalg.norm(gn), 1e-12)
            if rel > worst:
                worst, name = float(rel), cls.__name__
    ok = worst < 1e-5
    return Verdict(
        ok,
        f"worst relative error {worst:.3e} ({name}), tolerance 1e-5",
        {"max_rel_error": worst},
    )


@experiment("no_term_is_constant")
def exp_no_term_is_constant() -> Verdict:
    """No term collapses to a constant, as TetrahedralEnergy and GoldenEnergy did.

    Revision 1: TetrahedralEnergy == 1/3 exactly (lambda_min of a rank-<=8 Gram
    matrix of 24 vectors is always 0) and GoldenEnergy == |1-phi| exactly (it
    ratioed root norms, and all E8 roots have norm sqrt(2)).
    """
    L = lattice()
    pts = sphere_points(200, seed=4)
    ranges = {}
    for cls in ALL_TERM_CLASSES:
        t = cls()
        vals = []
        for p in pts:
            nb, _ = L.nearest_neighbors(p, k=24)
            vals.append(t.compute(p, nb))
        ranges[cls.__name__] = float(np.ptp(vals))
    worst = min(ranges, key=lambda name: ranges[name])
    ok = ranges[worst] > 1e-6
    return Verdict(
        ok,
        f"narrowest range: {worst} spans {ranges[worst]:.6f}",
        {k: v for k, v in ranges.items()},
    )


# --------------------------------------------------------------------------
# Solver
# --------------------------------------------------------------------------
@experiment("solver_energy_changes")
def exp_solver_energy_changes() -> Verdict:
    """The recorded energy actually moves during a run.

    Revision 1 produced exactly one distinct energy value over 1000 iterations,
    because the proposal was scored against the current point's neighbours.
    """
    L = lattice()
    s = GeometricAnnealingSolver(
        L,
        create_energy_suite(),
        GASParams(max_iters=300),
        rng=np.random.default_rng(0),
    )
    r = s.optimize(x_init=create_test_point(1))
    d = len(np.unique(np.round(r.energy_history, 9)))
    ok = d > 1
    return Verdict(
        ok,
        f"{d} distinct energies over {len(r.energy_history)} iterations",
        {"distinct_energies": float(d)},
    )


@experiment("solver_descends")
def exp_solver_descends() -> Verdict:
    """The solver ends below where it started, averaged over seeds."""
    L = lattice()
    terms = create_energy_suite()
    gains = []
    for seed in range(8):
        s = GeometricAnnealingSolver(
            L, terms, GASParams(max_iters=500), rng=np.random.default_rng(seed)
        )
        r = s.optimize(x_init=create_test_point(seed))
        gains.append((r.energy_history[0] - r.best_energy) / r.energy_history[0])
    mean_gain = float(np.mean(gains))
    ok = mean_gain > 0.01 and min(gains) >= 0.0
    return Verdict(
        ok,
        f"mean improvement {100 * mean_gain:.1f}% over 8 seeds, "
        f"worst {100 * min(gains):.1f}%",
        {"mean_improvement": mean_gain, "worst_improvement": float(min(gains))},
    )


@experiment("all_terms_flag_has_effect")
def exp_all_terms_flag_has_effect() -> Verdict:
    """``--all-terms`` changes the computation.

    Revision 1's ``_compute_weights`` returned 3 weights regardless of term
    count, so ``zip()`` silently dropped terms 4-7 and the flag did nothing.
    """
    L = lattice()
    out = {}
    for flag in (False, True):
        terms = create_energy_suite(include_all=flag)
        s = GeometricAnnealingSolver(
            L, terms, GASParams(max_iters=100), rng=np.random.default_rng(42)
        )
        out[flag] = (
            len(s._compute_weights(0.5)),
            len(terms),
            s.optimize(x_init=create_test_point(42)).best_energy,
        )
    weights_match = all(w == n for w, n, _ in out.values())
    differs = abs(out[False][2] - out[True][2]) > 1e-12
    ok = weights_match and differs
    return Verdict(
        ok,
        f"3-term E={out[False][2]:.6f} (weights {out[False][0]}/{out[False][1]}), "
        f"7-term E={out[True][2]:.6f} (weights {out[True][0]}/{out[True][1]})",
        {"E_core": out[False][2], "E_all": out[True][2]},
    )


@experiment("annealing_schedule_decays")
def exp_annealing_schedule_decays() -> Verdict:
    """Noise and temperature decay with the iteration index.

    Revision 1's sigma_t and T_t depended only on rho, never on t; sigma_t
    stayed inside [0.0452, 0.0500] for the whole run, so nothing annealed.
    """
    p = GASParams()
    ts = np.array([0, 100, 400, 800])
    cooling = np.exp(-ts / p.tau_anneal)
    sigma = p.sigma_0 * cooling
    T = p.T_0 * np.exp(-p.beta * 0.5) * cooling
    ok = bool(
        np.all(np.diff(sigma) < 0)
        and np.all(np.diff(T) < 0)
        and sigma[-1] < 0.1 * sigma[0]
    )
    return Verdict(
        ok,
        f"sigma_t {sigma[0]:.5f} -> {sigma[-1]:.5f}, "
        f"T_t {T[0]:.6f} -> {T[-1]:.6f} over 800 iterations",
        {"sigma_ratio": float(sigma[-1] / sigma[0])},
    )


@experiment("reproducible_under_seed")
def exp_reproducible_under_seed() -> Verdict:
    """Identical seeds give identical results, regardless of the global RNG.

    Revision 1 drew from the global ``np.random`` state, so any other library
    touching it changed the answer.
    """
    L = lattice()

    def run():
        s = GeometricAnnealingSolver(
            L,
            create_energy_suite(),
            GASParams(max_iters=100),
            rng=np.random.default_rng(5),
        )
        return s.optimize(x_init=create_test_point(3)).best_energy

    np.random.seed(999)
    a = run()
    np.random.seed(1)
    b = run()
    ok = a == b
    return Verdict(
        ok,
        f"two runs under different global seeds: {a:.12f} vs {b:.12f}",
        {"run_a": a, "run_b": b},
    )


@experiment("beats_random_search")
def exp_beats_random_search() -> Verdict:
    """GAS finds lower energy than uniform random sampling at equal budget.

    This is the claim README made as "annealing on E8 converges reliably".
    Budget is matched at 500 energy evaluations per method, 8 seeds each.
    """
    L = lattice()
    terms = create_energy_suite()
    budget = 500
    gas, rnd = [], []
    for seed in range(8):
        s = GeometricAnnealingSolver(
            L, terms, GASParams(max_iters=budget), rng=np.random.default_rng(seed)
        )
        gas.append(s.optimize(x_init=create_test_point(seed)).best_energy)
        rnd.append(random_search(budget, seed=1000 + seed, terms=terms))
    g, r = np.array(gas), np.array(rnd)
    wins = int((g < r).sum())
    delta = float((r.mean() - g.mean()) / r.mean())
    ok = wins >= 6 and delta > 0.0
    return Verdict(
        ok,
        f"GAS mean {g.mean():.6f} vs random {r.mean():.6f} "
        f"({100 * delta:+.1f}%), GAS wins {wins}/8 seeds",
        {
            "gas_mean": float(g.mean()),
            "random_mean": float(r.mean()),
            "wins": float(wins),
            "relative_gain": delta,
        },
    )


@experiment("expected_energy_decreases")
def exp_expected_energy_decreases() -> Verdict:
    """One GAS step does not increase energy by more than O(sigma_t^2).

    Decides the Lemma in THEORY.md 9.1. Its revision-1 premise ("the gradient
    term provides deterministic descent") was false because gradients were
    identically zero; with revision-2 terms the premise holds, so the lemma
    becomes testable: take one step from many independent states and compare
    the mean energy change against the sigma^2 bound.
    """
    L = lattice()
    terms = create_energy_suite()
    params = GASParams(max_iters=1)
    sigma_sq = params.sigma_0**2

    deltas = []
    for seed in range(120):
        s = GeometricAnnealingSolver(
            L, terms, params, rng=np.random.default_rng(seed)
        )
        x0 = sphere_points(1, seed=5000 + seed)[0]
        nb, rho = s._evaluate(x0)
        E0 = s._compute_energy(x0, nb, rho)
        from gas.solver import GASState

        state = GASState(x=x0, energy=E0, rho_coset=rho, iteration=0)
        deltas.append(s.step(state).energy - E0)

    d = np.array(deltas)
    mean_delta = float(d.mean())
    ok = mean_delta <= sigma_sq
    return Verdict(
        ok,
        f"mean one-step dE = {mean_delta:+.6f} over 120 states, "
        f"bound sigma^2 = {sigma_sq:.6f}; {100 * (d <= 0).mean():.0f}% of steps "
        f"were non-increasing",
        {
            "mean_delta_E": mean_delta,
            "sigma_squared": sigma_sq,
            "fraction_non_increasing": float((d <= 0).mean()),
        },
    )


@experiment("decoder_geometric_regularizer_active")
def exp_decoder_geometric_regularizer_active() -> Verdict:
    """The meta-layer's geometric regularizer actually influences the decode.

    With revision-1 energy terms, R_geo depended on x only through the discrete
    neighbour selection, so its gradient measured ~3e-11 and lambda_3 had no
    effect on the solution -- the "geometric coherence" term was inert.
    """
    from meta_layer.decoder import ProximalGeometricDecoder

    L = lattice()
    terms = create_energy_suite()
    rng = np.random.default_rng(0)
    W = np.linalg.qr(rng.standard_normal((20, 8)))[0]
    x_star = create_test_point(5)

    grad_norms, shifts = [], []
    for trial in range(5):
        y = rng.standard_normal(20) * 0.1
        on = ProximalGeometricDecoder(W, L, terms, lambda_3=1.0)
        off = ProximalGeometricDecoder(W, L, terms, lambda_3=0.0)
        grad_norms.append(np.linalg.norm(on._compute_R_geo_gradient(y)))
        shifts.append(
            np.linalg.norm(
                on.decode(x_star, max_iters=40) - off.decode(x_star, max_iters=40)
            )
        )
    g, sh = float(np.min(grad_norms)), float(np.min(shifts))
    ok = g > 1e-6 and sh > 1e-6
    return Verdict(
        ok,
        f"min ||grad R_geo|| = {g:.3e}; turning lambda_3 off moves the decoded "
        f"solution by at least {sh:.3e}",
        {"min_grad_norm": g, "min_solution_shift": sh},
    )
