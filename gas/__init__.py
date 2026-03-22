"""
Geometric Annealing Solver (GAS) for E8 Optimization

Core module providing the E8 lattice structure, geometric energy terms,
and the annealing-based solver.
"""

from .lattice import E8Lattice, create_test_point
from .energy_terms import (
    EnergyTerm,
    OctahedralEnergy,
    TetrahedralEnergy,
    GoldenEnergy,
    SquareEnergy,
    HexagonalEnergy,
    DodecahedralEnergy,
    create_energy_suite,
)
from .solver import GeometricAnnealingSolver, GASParams, GASState

__all__ = [
    "E8Lattice",
    "create_test_point",
    "EnergyTerm",
    "OctahedralEnergy",
    "TetrahedralEnergy",
    "GoldenEnergy",
    "SquareEnergy",
    "HexagonalEnergy",
    "DodecahedralEnergy",
    "create_energy_suite",
    "GeometricAnnealingSolver",
    "GASParams",
    "GASState",
]
