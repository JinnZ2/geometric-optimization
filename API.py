###examples!!###

"""
End-to-end G-Opt workflow example
"""
from gas import E8Lattice, GeometricAnnealingSolver, GASParams
from gas.energy_terms import OctahedralEnergy, TetrahedralEnergy, GoldenEnergy
from meta_layer import WPCAEncoder, ProximalGeometricDecoder

# 1. Initialize E₈ lattice
lattice = E8Lattice()

# 2. Define energy terms
energy_terms = [
    OctahedralEnergy(),
    TetrahedralEnergy(),
    GoldenEnergy()
]

# 3. Configure solver
params = GASParams(
    k_neighbors=24,
    max_iters=1000,
    rho_min=0.6
)

# 4. Create solver
solver = GeometricAnnealingSolver(lattice, energy_terms, params)

# 5. Optimize in E₈ space
result = solver.optimize()

print(f"Converged: {result.converged}")
print(f"Final energy: {result.energy:.6f}")
print(f"Coset density: {result.rho_coset:.3f}")

# 6. Decode back to N-dimensional space
# (Assuming we have W from prior W-PCA encoding)
decoder = ProximalGeometricDecoder(W, lattice, energy_terms)
y_optimal = decoder.decode(result.x)

print(f"Optimal solution (N-space): {y_optimal}")
