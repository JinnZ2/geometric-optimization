"""
Geometric Annealing Solver (GAS) for E8 Optimization

This package implements a novel optimization framework that replaces
linear cost minimization with geometric coherence maximization on the
E8 exceptional Lie group lattice.

Key Components:
- E8Lattice: 240-root system with D8/coset decomposition
- Energy Terms: Spectral formulations (octahedral, tetrahedral, golden)
- GAS Solver: Adaptive annealing with phi-folding transformation

Example:
>>> from gas import E8Lattice, GeometricAnnealingSolver, GASParams
>>> from gas import create_energy_suite
>>>
>>> lattice = E8Lattice()
>>> energy_terms = create_energy_suite()
>>> params = GASParams(max_iters=100)
>>> solver = GeometricAnnealingSolver(lattice, energy_terms, params)
>>>
>>> result = solver.optimize()
>>> print(f"Final energy: {result.energy:.6f}")
>>> print(f"phi-alignment: {result.rho_coset:.3f}")

For detailed mathematical theory, see THEORY.md in the repository root.

References:
- Viazovska (2016): E8 optimal sphere packing
- Conway & Sloane (1988): Sphere Packings, Lattices and Groups
- Kirkpatrick et al. (1983): Simulated Annealing
"""

from gas.lattice import E8Lattice, create_test_point
from gas.energy_terms import (
    EnergyTerm,
    OctahedralEnergy,
    TetrahedralEnergy,
    GoldenEnergy,
    SquareEnergy,
    HexagonalEnergy,
    DodecahedralEnergy,
    create_energy_suite,
)
from gas.solver import (
    GeometricAnnealingSolver,
    GASParams,
    GASState,
)

__version__ = "0.1.0"
__author__ = "JinnZ2 Collaborative Research"
__license__ = "Apache 2.0"

__all__ = [
    # Lattice
    "E8Lattice",
    "create_test_point",
    # Energy Terms
    "EnergyTerm",
    "OctahedralEnergy",
    "TetrahedralEnergy",
    "GoldenEnergy",
    "SquareEnergy",
    "HexagonalEnergy",
    "DodecahedralEnergy",
    "create_energy_suite",
    # Solver
    "GeometricAnnealingSolver",
    "GASParams",
    "GASState",
]
