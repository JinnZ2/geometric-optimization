"""
E₈ Lattice Structure with D₈/Coset Decomposition
"""
import numpy as np
from scipy.spatial import KDTree
from dataclasses import dataclass

@dataclass
class E8Lattice:
    """Complete E₈ root system (240 roots, norm²=2)"""
    
    def __init__(self):
        self.d8_roots = self._generate_d8_roots()      # 112 Cartesian
        self.coset_roots = self._generate_coset()      # 128 half-integer
        self.all_roots = np.vstack([self.d8_roots, self.coset_roots])
        self.kdtree = KDTree(self.all_roots)
        
        # Precompute coset membership mask
        self.is_coset = np.zeros(240, dtype=bool)
        self.is_coset[112:] = True
    
    def _generate_d8_roots(self) -> np.ndarray:
        """Generate 112 D₈ roots: {±eᵢ ± eⱼ | i<j}"""
        roots = []
        # Type 1: ±eᵢ ± eⱼ (i ≠ j) - 56 pairs × 2 signs = 112 roots
        for i in range(8):
            for j in range(i+1, 8):
                for sign_i in [1, -1]:
                    for sign_j in [1, -1]:
                        root = np.zeros(8)
                        root[i] = sign_i
                        root[j] = sign_j
                        roots.append(root)
        return np.array(roots)
    
    def _generate_coset(self) -> np.ndarray:
        """Generate 128 half-integer coset roots: (±½)⁸ with even # of -½"""
        roots = []
        for bits in range(256):  # All 8-bit patterns
            root = np.array([1 if bits & (1 << i) else -1 
                           for i in range(8)]) / 2
            # Keep only even parity (even number of -½ components)
            if np.sum(root < 0) % 2 == 0:
                roots.append(root)
        return np.array(roots)
    
    def nearest_neighbors(self, x: np.ndarray, k: int = 24):
        """Find k-nearest E₈ roots to point x"""
        distances, indices = self.kdtree.query(x, k=k)
        return self.all_roots[indices], indices
    
    def coset_density(self, indices: np.ndarray) -> float:
        """Calculate ρ_coset for neighborhood"""
        return np.mean(self.is_coset[indices])
