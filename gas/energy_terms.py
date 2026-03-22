"""
Geometric Energy Terms (Spectral Formulations)
"""
import numpy as np
from abc import ABC, abstractmethod

class EnergyTerm(ABC):
    """Base class for geometric energy terms"""
    
    @abstractmethod
    def compute(self, x: np.ndarray, neighbors: np.ndarray) -> float:
        """Calculate energy for state x given neighborhood"""
        pass
    
    def gradient(self, x: np.ndarray, neighbors: np.ndarray, 
                 eps: float = 1e-3) -> np.ndarray:
        """Finite difference gradient approximation"""
        grad = np.zeros(8)
        E0 = self.compute(x, neighbors)
        for i in range(8):
            x_plus = x.copy()
            x_plus[i] += eps
            x_minus = x.copy()
            x_minus[i] -= eps
            grad[i] = (self.compute(x_plus, neighbors) - 
                      self.compute(x_minus, neighbors)) / (2 * eps)
        return grad


class OctahedralEnergy(EnergyTerm):
    """Spectral octahedral term: λ₁(G - I)"""
    
    def compute(self, x: np.ndarray, neighbors: np.ndarray) -> float:
        # Normalize neighbors to unit vectors
        norms = np.linalg.norm(neighbors, axis=1, keepdims=True)
        unit_neighbors = neighbors / (norms + 1e-10)
        
        # Gram matrix
        G = unit_neighbors @ unit_neighbors.T
        
        # Deviation from identity (perfect orthogonality)
        deviation = G - np.eye(len(G))
        eigenvalues = np.linalg.eigvalsh(deviation)
        
        return np.max(np.abs(eigenvalues))


class TetrahedralEnergy(EnergyTerm):
    """Spectral tetrahedral term: |λ_min(G) + 1/3|"""
    
    def compute(self, x: np.ndarray, neighbors: np.ndarray) -> float:
        norms = np.linalg.norm(neighbors, axis=1, keepdims=True)
        unit_neighbors = neighbors / (norms + 1e-10)
        
        G = unit_neighbors @ unit_neighbors.T
        eigenvalues = np.linalg.eigvalsh(G)
        
        # Tetrahedral angle cosine: -1/3
        return np.abs(np.min(eigenvalues) + 1/3)


class GoldenEnergy(EnergyTerm):
    """φ-ratio alignment energy"""
    
    def compute(self, x: np.ndarray, neighbors: np.ndarray) -> float:
        phi = (1 + np.sqrt(5)) / 2
        
        norms = np.linalg.norm(neighbors, axis=1)
        if len(norms) < 2:
            return 0.0
        
        # Pairwise ratios
        ratios = norms[:, None] / (norms[None, :] + 1e-10)
        
        # Deviation from φ
        deviation = np.abs(ratios - phi)
        
        # Use median to be robust to outliers
        return np.median(deviation)


class SquareEnergy(EnergyTerm):
    """Square symmetry energy: deviation from 90-degree angles"""

    def compute(self, x: np.ndarray, neighbors: np.ndarray) -> float:
        norms = np.linalg.norm(neighbors, axis=1, keepdims=True)
        unit_neighbors = neighbors / (norms + 1e-10)
        dots = unit_neighbors @ unit_neighbors.T
        # Ideal square angles: 0 (90 degrees) between distinct pairs
        np.fill_diagonal(dots, 0.0)
        return np.mean(dots ** 2)


class HexagonalEnergy(EnergyTerm):
    """Hexagonal symmetry energy: deviation from 60-degree packing"""

    def compute(self, x: np.ndarray, neighbors: np.ndarray) -> float:
        norms = np.linalg.norm(neighbors, axis=1, keepdims=True)
        unit_neighbors = neighbors / (norms + 1e-10)
        dots = unit_neighbors @ unit_neighbors.T
        np.fill_diagonal(dots, 0.0)
        # Hexagonal packing angle cos(60°) = 0.5
        return np.mean((np.abs(dots) - 0.5) ** 2)


class DodecahedralEnergy(EnergyTerm):
    """Dodecahedral symmetry energy: deviation from icosahedral angles"""

    def compute(self, x: np.ndarray, neighbors: np.ndarray) -> float:
        phi = (1 + np.sqrt(5)) / 2
        # Icosahedral angle cosine: 1/sqrt(5) ≈ 1/φ²
        target = 1.0 / (phi ** 2)
        norms = np.linalg.norm(neighbors, axis=1, keepdims=True)
        unit_neighbors = neighbors / (norms + 1e-10)
        dots = unit_neighbors @ unit_neighbors.T
        np.fill_diagonal(dots, 0.0)
        return np.mean((np.abs(dots) - target) ** 2)


def create_energy_suite(include_all: bool = False):
    """Create the standard set of geometric energy terms.

    Args:
        include_all: If True, include square, hexagonal, and dodecahedral
            terms in addition to the core three.

    Returns:
        List of EnergyTerm instances.
    """
    terms = [OctahedralEnergy(), TetrahedralEnergy(), GoldenEnergy()]
    if include_all:
        terms.extend([SquareEnergy(), HexagonalEnergy(), DodecahedralEnergy()])
    return terms
