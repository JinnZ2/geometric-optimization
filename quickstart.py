#!/usr/bin/env python3
"""
Quick Start Example for Geometric Optimization Framework

This script demonstrates basic usage of the GAS solver to find
optimal configurations on the E8 lattice.

Run with: python quickstart.py
"""

import numpy as np
import sys

sys.path.insert(0, ".")

from gas.lattice import E8Lattice, create_test_point
from gas.energy_terms import create_energy_suite
from gas.solver import GeometricAnnealingSolver, GASParams


def main():
    """Run basic optimization example."""
    print("=" * 70)
    print(" Geometric Optimization Framework - Quick Start")
    print("=" * 70)
    print()

    # 1. Initialize E8 lattice
    print("1. Initializing E8 lattice...")
    lattice = E8Lattice()
    print(f"   Loaded {len(lattice.all_roots)} roots ({len(lattice.d8_roots)} D8 + {len(lattice.coset_roots)} coset)")
    print()

    # 2. Create energy terms
    print("2. Creating geometric energy terms...")
    energy_terms = create_energy_suite(include_all=False)
    print(f"   Loaded {len(energy_terms)} energy terms:")
    for term in energy_terms:
        print(f"     - {term.__class__.__name__}")
    print()

    # 3. Configure solver
    print("3. Configuring GAS solver...")
    params = GASParams(
        max_iters=300,
        rho_min=0.5,
        tau_E=1e-4,
    )
    print(f"   Max iterations: {params.max_iters}")
    print(f"   Target rho_coset: {params.rho_min}")
    print()

    # 4. Initialize solver
    print("4. Initializing solver...")
    solver = GeometricAnnealingSolver(lattice, energy_terms, params)
    print("   phi-rotation matrix constructed")
    print("   Adaptive annealing schedules configured")
    print()

    # 5. Run optimization
    print("5. Running geometric annealing...")
    print()

    # Reproducible initial state
    x_init = create_test_point(seed=42)

    # Progress callback
    def progress_callback(state):
        if state.iteration % 50 == 0:
            print(
                f"   Iteration {state.iteration:3d}: "
                f"E={state.energy:.6f}, "
                f"rho={state.rho_coset:.3f}"
            )

    result = solver.optimize(x_init=x_init, callback=progress_callback)

    print()
    if result.converged:
        print("   Converged!")
    else:
        print("   Maximum iterations reached")
    print()

    # 6. Report results
    print("=" * 70)
    print(" Results")
    print("=" * 70)
    print()
    print(f"Final Energy:       {result.energy:.6f}")
    print(f"Coset Density (rho): {result.rho_coset:.3f}")
    print(f"Iterations:         {result.iteration}")
    print()
    print(f"Energy Reduction:   {(1 - result.energy / result.energy_history[0]) * 100:.1f}%")
    print(f"Initial Energy:     {result.energy_history[0]:.6f}")
    print(f"Minimum Energy:     {min(result.energy_history):.6f}")
    print(f"Final Energy:       {result.energy_history[-1]:.6f}")
    print()

    # 7. Geometric interpretation
    print("=" * 70)
    print(" Geometric Interpretation")
    print("=" * 70)
    print()

    neighbors, indices = lattice.nearest_neighbors(result.x, k=24)
    n_d8 = np.sum(~lattice.is_coset[indices])
    n_coset = np.sum(lattice.is_coset[indices])

    print(f"Final neighborhood composition:")
    print(f"  D8 roots (Cartesian):      {n_d8}/24 ({n_d8 / 24 * 100:.0f}%)")
    print(f"  Coset roots (phi-rich):    {n_coset}/24 ({n_coset / 24 * 100:.0f}%)")
    print()

    if result.rho_coset > 0.6:
        print("  -> System achieved EXCEPTIONAL PHASE")
        print("  -> Configuration exhibits phi-rich quasi-crystalline order")
        print("  -> High geometric information density")
    elif result.rho_coset < 0.4:
        print("  -> System in CARTESIAN PHASE")
        print("  -> Configuration exhibits orthogonal linear order")
        print("  -> Lower geometric information density")
    else:
        print("  -> System in TRANSITION REGION")
        print("  -> Mixed geometric order")
        print("  -> Balanced information density")
    print()

    # 8. Conclusion
    print("=" * 70)
    print(" Quick Start Complete")
    print("=" * 70)
    print()
    print("Next steps:")
    print("  - Read THEORY.md for mathematical foundations")
    print("  - Read GUIDE.md for the E8 <-> Rosetta connection")
    print("  - See bridges/rosetta-fieldlink.json for the polyhedral mapping")
    print()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\nInterrupted by user")
        sys.exit(0)
    except Exception as e:
        print(f"\n\nError: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)
