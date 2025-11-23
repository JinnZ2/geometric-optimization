“””
Geometric Annealing Solver (GAS) for E₈ Optimization

This package implements a novel optimization framework that replaces
linear cost minimization with geometric coherence maximization on the
E₈ exceptional Lie group lattice.

Key Components:
- E8Lattice: 240-root system with D₈/coset decomposition
- Energy Terms: Spectral formulations (octahedral, tetrahedral, golden)
- GAS Solver: Adaptive annealing with φ-folding transformation

Example:
>>> from gas import E8Lattice, GeometricAnnealingSolver, GASParams
>>> from gas import create_energy_suite
>>>
>>> lattice = E8Lattice()
>>> energy_terms = create_energy_suite()
>>> solver = GeometricAnnealingSolver(lattice, energy_terms)
>>>
>>> result = solver.optimize()
>>> print(f”Final energy: {result.energy:.6f}”)
>>> print(f”φ-alignment: {result.rho_coset:.3f}”)

For detailed mathematical theory, see THEORY.md in the repository root.

References:
- Viazovska (2016): E₈ optimal sphere packing
- Conway & Sloane (1988): Sphere Packings, Lattices and Groups
- Kirkpatrick et al. (1983): Simulated Annealing
“””

from .lattice import E8Lattice, create_test_point
from .energy_terms import (
EnergyTerm,
OctahedralEnergy,
TetrahedralEnergy,
GoldenEnergy,
SquareEnergy,
HexagonalEnergy,
DodecahedralEnergy,
create_energy_suite
)
from .solver import (
GeometricAnnealingSolver,
GASParams,
GASState
)

**version** = “0.1.0”
**author** = “JinnZ2 Collaborative Research”
**license** = “Apache 2.0”

**all** = [
# Lattice
“E8Lattice”,
“create_test_point”,

```
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
```

]
