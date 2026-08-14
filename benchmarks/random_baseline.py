"""
GAS against the trivial baseline: uniform random sampling of the sphere.

An optimizer that cannot beat random sampling at an equal evaluation budget has
not been shown to optimize anything. This script produces the measurement that
decides claim C11 in ``validation/claims.json``.

Run::

    python -m benchmarks.random_baseline --budget 500 --seeds 8
"""

from __future__ import annotations

import argparse
import numpy as np

from gas.lattice import E8Lattice, create_test_point
from gas.energy_terms import create_energy_suite
from gas.solver import GeometricAnnealingSolver, GASParams


def random_search(lattice, probe, budget: int, seed: int) -> float:
    """Best energy found by drawing ``budget`` uniform points on the sphere."""
    rng = np.random.default_rng(seed)
    best = np.inf
    for _ in range(budget):
        p = rng.standard_normal(8)
        p = p / np.linalg.norm(p) * np.sqrt(2)
        nb, idx = lattice.nearest_neighbors(p, k=24)
        best = min(best, probe._compute_energy(p, nb, lattice.coset_density(idx)))
    return float(best)


def multistart_gas(lattice, terms, budget: int, seed: int, restarts: int) -> float:
    """Best energy from ``restarts`` short GAS runs sharing the same budget."""
    rng = np.random.default_rng(seed)
    per = max(1, budget // restarts)
    best = np.inf
    for _ in range(restarts):
        s = GeometricAnnealingSolver(
            lattice,
            terms,
            GASParams(max_iters=per, tau_anneal=max(per / 2, 1.0)),
            rng=rng,
        )
        best = min(best, s.optimize().best_energy)
    return float(best)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="random_baseline")
    ap.add_argument("--budget", type=int, default=500, help="energy evals per method")
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--restarts", type=int, default=16, help="multi-start GAS restarts")
    args = ap.parse_args(argv)

    lattice = E8Lattice()
    terms = create_energy_suite()
    probe = GeometricAnnealingSolver(
        lattice, terms, GASParams(), rng=np.random.default_rng(0)
    )

    gas, multi, rand = [], [], []
    for seed in range(args.seeds):
        s = GeometricAnnealingSolver(
            lattice, terms, GASParams(max_iters=args.budget),
            rng=np.random.default_rng(seed),
        )
        gas.append(s.optimize(x_init=create_test_point(seed)).best_energy)
        multi.append(
            multistart_gas(lattice, terms, args.budget, seed, args.restarts)
        )
        rand.append(random_search(lattice, probe, args.budget, 1000 + seed))

    g, m, r = np.array(gas), np.array(multi), np.array(rand)
    print(f"budget = {args.budget} energy evaluations, {args.seeds} seeds\n")
    print(f"{'method':<22}{'mean':>12}{'std':>12}{'wins vs random':>16}")
    print("-" * 62)
    print(f"{'random sampling':<22}{r.mean():>12.6f}{r.std():>12.6f}{'--':>16}")
    print(
        f"{'GAS (single start)':<22}{g.mean():>12.6f}{g.std():>12.6f}"
        f"{f'{int((g < r).sum())}/{args.seeds}':>16}"
    )
    print(
        f"{'GAS (multi start)':<22}{m.mean():>12.6f}{m.std():>12.6f}"
        f"{f'{int((m < r).sum())}/{args.seeds}':>16}"
    )
    single = 100 * (r.mean() - g.mean()) / r.mean()
    multi_pct = 100 * (r.mean() - m.mean()) / r.mean()
    print(
        f"\nrelative to random: single start {single:+.1f}%"
        f", multi start {multi_pct:+.1f}%"
    )
    print(
        "\nA method must win a clear majority of seeds to count as beating the\n"
        "baseline. See claim C11 in VALIDATION.md."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
