"""
End-to-end G-Opt workflow example

Demonstrates the full pipeline:
  1. Initialize E8 lattice
  2. Optimize in 8D (GAS solver)
  3. Decode back to N-dimensional space (ProximalGeometricDecoder)
"""

import numpy as np
from gas import E8Lattice, GeometricAnnealingSolver, GASParams
from gas.energy_terms import OctahedralEnergy, TetrahedralEnergy, GoldenEnergy
from meta_layer import ProximalGeometricDecoder

# 1. Initialize E8 lattice
lattice = E8Lattice()

# 2. Define energy terms
energy_terms = [
    OctahedralEnergy(),
    TetrahedralEnergy(),
    GoldenEnergy(),
]

# 3. Configure solver
params = GASParams(
    k_neighbors=24,
    max_iters=1000,
    rho_min=0.6,
)

# 4. Create solver and optimize in E8 space
solver = GeometricAnnealingSolver(lattice, energy_terms, params)
result = solver.optimize()

print(f"Converged: {result.converged}")
print(f"Final energy: {result.energy:.6f}")
print(f"Coset density: {result.rho_coset:.3f}")

# 5. Decode back to N-dimensional space
# W is an (N, 8) projection matrix that maps N-dimensional problems into E8.
# In a real application, W comes from domain-specific dimensionality reduction
# (e.g., PCA on your problem's feature space projected to 8 components).
# Here we demonstrate with a random 20-dimensional problem:
N = 20
rng = np.random.RandomState(42)
W = rng.randn(N, 8)
W = np.linalg.qr(W.T)[0].T  # Orthonormalize columns for stability

decoder = ProximalGeometricDecoder(W, lattice, energy_terms)
y_optimal = decoder.decode(result.x)

print(f"Optimal solution (N={N} space): shape={y_optimal.shape}")
print(f"  L2 norm: {np.linalg.norm(y_optimal):.4f}")
print(f"  Sparsity (|y| < 0.01): {np.sum(np.abs(y_optimal) < 0.01)}/{N}")
