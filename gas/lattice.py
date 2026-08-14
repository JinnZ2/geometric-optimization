"""
E8 Lattice Structure with D8/Coset Decomposition
"""

import numpy as np
from scipy.spatial import KDTree
from dataclasses import dataclass


@dataclass
class E8Lattice:
    """Complete E8 root system (240 roots, norm^2=2)"""

    def __init__(self):
        self.d8_roots = self._generate_d8_roots()  # 112 Cartesian
        self.coset_roots = self._generate_coset()  # 128 half-integer
        self.all_roots = np.vstack([self.d8_roots, self.coset_roots])
        self.kdtree = KDTree(self.all_roots)

        # Precompute coset membership mask
        self.is_coset = np.zeros(240, dtype=bool)
        self.is_coset[112:] = True

    def _generate_d8_roots(self) -> np.ndarray:
        """Generate 112 D8 roots: {+/-ei +/- ej | i<j}"""
        roots = []
        for i in range(8):
            for j in range(i + 1, 8):
                for sign_i in [1, -1]:
                    for sign_j in [1, -1]:
                        root = np.zeros(8)
                        root[i] = sign_i
                        root[j] = sign_j
                        roots.append(root)
        return np.array(roots)

    def _generate_coset(self) -> np.ndarray:
        """Generate 128 half-integer coset roots: (+/-1/2)^8 with even # of -1/2"""
        roots = []
        for bits in range(256):  # All 8-bit patterns
            root = np.array([1 if bits & (1 << i) else -1 for i in range(8)]) / 2
            # Keep only even parity (even number of -1/2 components)
            if np.sum(root < 0) % 2 == 0:
                roots.append(root)
        return np.array(roots)

    def nearest_neighbors(self, x: np.ndarray, k: int = 24):
        """Find k-nearest E8 roots to point x"""
        distances, indices = self.kdtree.query(x, k=k)
        return self.all_roots[indices], indices

    def coset_density(self, indices: np.ndarray) -> float:
        """Calculate rho_coset for neighborhood"""
        return float(np.mean(self.is_coset[indices]))


def create_test_point(seed: int = 42) -> np.ndarray:
    """Create a reproducible random test point on the E8 norm-sqrt(2) sphere."""
    rng = np.random.RandomState(seed)
    x = rng.randn(8)
    point: np.ndarray = x / np.linalg.norm(x) * np.sqrt(2)
    return point
