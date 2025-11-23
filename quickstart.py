#!/usr/bin/env python3
“””
Quick Start Example for Geometric Optimization Framework

This script demonstrates basic usage of the GAS solver to find
optimal configurations on the E₈ lattice.

Run with: python quickstart.py
“””

import numpy as np
import sys

# Note: When installed as package, imports will be:

# from gas import E8Lattice, GeometricAnnealingSolver, GASParams

# from gas import create_energy_suite

# For standalone execution, adjust path

sys.path.insert(0, ‘.’)

try:
from gas_lattice import E8Lattice
from gas_energy_terms import create_energy_suite
from gas_solver import GeometricAnnealingSolver, GASParams
except ImportError:
print(“Error: Could not import GAS modules”)
print(“Make sure gas_lattice.py, gas_energy_terms.py, and gas_solver.py are in the current directory”)
sys.exit(1)

def main():
“”“Run basic optimization example.”””
print(”=”*70)
print(” Geometric Optimization Framework - Quick Start”)
print(”=”*70)
print()

```
# 1. Initialize E₈ lattice
print("1. Initializing E₈ lattice...")
lattice = E8Lattice()
print(f"   ✓ {lattice}")
print()

# 2. Create energy terms
print("2. Creating geometric energy terms...")
energy_terms = create_energy_suite(include_all=False)
print(f"   ✓ Loaded {len(energy_terms)} energy terms:")
for term in energy_terms:
    print(f"     - {term}")
print()

# 3. Configure solver
print("3. Configuring GAS solver...")
params = GASParams(
    max_iters=300,
    rho_min=0.5,
    tau_E=1e-4,
    verbose=False  # Will handle output manually
)
print(f"   ✓ Max iterations: {params.max_iters}")
print(f"   ✓ Target ρ_coset: {params.rho_min}")
print()

# 4. Initialize solver
print("4. Initializing solver...")
solver = GeometricAnnealingSolver(lattice, energy_terms, params)
print("   ✓ φ-rotation matrix constructed")
print("   ✓ Adaptive annealing schedules configured")
print()

# 5. Run optimization
print("5. Running geometric annealing...")
print()

# Random initial state
np.random.seed(42)  # Reproducibility
x_init = np.random.randn(8)
x_init = x_init / np.linalg.norm(x_init) * np.sqrt(2)

# Progress callback
def progress_callback(state):
    if state.iteration % 50 == 0:
        print(f"   Iteration {state.iteration:3d}: "
              f"E={state.energy:.6f}, "
              f"ρ={state.rho_coset:.3f}, "
              f"accept={state.acceptance_rate:.1%}")

# Optimize
result = solver.optimize(x_init=x_init, callback=progress_callback, verbose=False)

print()
if result.converged:
    print("   ✓ Converged!")
else:
    print("   ⚠ Maximum iterations reached")
print()

# 6. Report results
print("="*70)
print(" Results")
print("="*70)
print()
print(f"Final Energy:       {result.energy:.6f}")
print(f"Coset Density (ρ):  {result.rho_coset:.3f}")
print(f"Iterations:         {result.iteration}")
print(f"Acceptance Rate:    {result.acceptance_rate:.1%}")
print()
print(f"Energy Reduction:   {(1 - result.energy/result.energy_history[0])*100:.1f}%")
print(f"Initial Energy:     {result.energy_history[0]:.6f}")
print(f"Minimum Energy:     {min(result.energy_history):.6f}")
print(f"Final Energy:       {result.energy_history[-1]:.6f}")
print()

# 7. Geometric interpretation
print("="*70)
print(" Geometric Interpretation")
print("="*70)
print()

# Get final neighborhood
neighbors, indices = lattice.nearest_neighbors(result.x, k=24)

# Count D₈ vs coset
n_d8 = np.sum(~lattice.is_coset[indices])
n_coset = np.sum(lattice.is_coset[indices])

print(f"Final neighborhood composition:")
print(f"  D₈ roots (Cartesian):      {n_d8}/24 ({n_d8/24*100:.0f}%)")
print(f"  Coset roots (φ-rich):      {n_coset}/24 ({n_coset/24*100:.0f}%)")
print()

if result.rho_coset > 0.6:
    print("  → System achieved EXCEPTIONAL PHASE")
    print("  → Configuration exhibits φ-rich quasi-crystalline order")
    print("  → High geometric information density")
elif result.rho_coset < 0.4:
    print("  → System in CARTESIAN PHASE")
    print("  → Configuration exhibits orthogonal linear order")
    print("  → Lower geometric information density")
else:
    print("  → System in TRANSITION REGION")
    print("  → Mixed geometric order")
    print("  → Balanced information density")
print()

# 8. Conclusion
print("="*70)
print(" Quick Start Complete")
print("="*70)
print()
print("Next steps:")
print("  • Explore examples/ for more detailed demonstrations")
print("  • Read THEORY.md for mathematical foundations")
print("  • See CONTRIBUTING.md to extend the framework")
print()
print("Framework ready for geometric optimization!")
print()
```

if **name** == “**main**”:
try:
main()
except KeyboardInterrupt:
print(”\n\nInterrupted by user”)
sys.exit(0)
except Exception as e:
print(f”\n\nError: {e}”)
import traceback
traceback.print_exc()
sys.exit(1)
