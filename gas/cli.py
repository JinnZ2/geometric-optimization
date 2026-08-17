"""
Command-line interface for the Geometric Annealing Solver.
"""

import argparse
import sys
import numpy as np

from .lattice import E8Lattice, create_test_point
from .energy_terms import create_energy_suite
from .solver import GeometricAnnealingSolver, GASParams


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="gopt-optimize",
        description=(
            "Geometric Annealing Solver (GAS) — "
            "optimize on the E8 lattice"
        ),
    )
    parser.add_argument(
        "--max-iters", type=int, default=1000,
        help="Maximum iterations (default: 1000)",
    )
    parser.add_argument(
        "--seed", type=int, default=None,
        help="Random seed for reproducibility",
    )
    parser.add_argument(
        "--rho-min", type=float, default=0.6,
        help="Minimum coset density target (default: 0.6)",
    )
    parser.add_argument(
        "--all-terms", action="store_true",
        help="Include all energy terms (square, hexagonal, "
             "dodecahedral, icosahedral)",
    )
    parser.add_argument(
        "--quiet", action="store_true",
        help="Suppress progress output",
    )
    args = parser.parse_args(argv)

    # The solver draws from an injected Generator, so seeding must go through
    # it -- setting the global np.random state would no longer reach the solver.
    rng = np.random.default_rng(args.seed)

    lattice = E8Lattice()
    energy_terms = create_energy_suite(include_all=args.all_terms)
    params = GASParams(max_iters=args.max_iters, rho_min=args.rho_min)
    solver = GeometricAnnealingSolver(lattice, energy_terms, params, rng=rng)

    x_init = create_test_point(
        args.seed if args.seed is not None else 42
    )

    def callback(state):
        if not args.quiet and state.iteration % 100 == 0:
            print(
                f"  iter {state.iteration:4d}  "
                f"E={state.energy:.6f}  "
                f"rho={state.rho_coset:.3f}"
            )

    if not args.quiet:
        print("Running GAS optimization...")

    result = solver.optimize(x_init=x_init, callback=callback)

    if not args.quiet:
        print()

    status = "CONVERGED" if result.converged else "MAX_ITERS"
    print(f"Status:      {status}")
    print(f"Iterations:  {result.iteration}")
    print(f"Energy:      {result.energy:.6f}")
    print(f"Rho_coset:   {result.rho_coset:.3f}")

    return 0 if result.converged else 1


if __name__ == "__main__":
    sys.exit(main())
