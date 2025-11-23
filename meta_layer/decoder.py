"""
Proximal-Geometric Descent (PGD) Inverse Solver
Decodes E₈ state back to N-dimensional actionable space
"""
import numpy as np
from scipy.optimize import minimize

class ProximalGeometricDecoder:
    """8→N inverse mapping with geometric regularization"""
    
    def __init__(self, 
                 W: np.ndarray,
                 lattice: 'E8Lattice',
                 energy_terms: List[EnergyTerm],
                 lambda_1: float = 0.01,  # L1 sparsity
                 lambda_2: float = 0.1,   # Cost model
                 lambda_3: float = 1.0):  # Geometric coherence
        self.W = W  # (N, 8) projection matrix
        self.lattice = lattice
        self.energy_terms = energy_terms
        self.lambda_1 = lambda_1
        self.lambda_2 = lambda_2
        self.lambda_3 = lambda_3
    
    def decode(self, x_star: np.ndarray, 
               cost_model: Optional[Callable] = None,
               max_iters: int = 500) -> np.ndarray:
        """
        Solve: min ||W^T y - x*||² + λ₁||y||₁ + λ₂·Cost(y) + λ₃·R_geo(y)
        """
        N = self.W.shape[0]
        
        # Warm start: least squares solution
        WtW = self.W.T @ self.W
        y_init = self.W @ np.linalg.solve(WtW + 1e-6*np.eye(8), x_star)
        
        def objective(y):
            # Fidelity term
            x_proj = self.W.T @ y
            fidelity = np.sum((x_proj - x_star)**2)
            
            # Sparsity term
            sparsity = self.lambda_1 * np.sum(np.abs(y))
            
            # Cost model term
            cost = 0.0
            if cost_model is not None:
                cost = self.lambda_2 * cost_model(y)
            
            # Geometric regularization
            geo_reg = self.lambda_3 * self._compute_R_geo(y)
            
            return fidelity + sparsity + cost + geo_reg
        
        def gradient(y):
            grad = np.zeros(N)
            
            # Fidelity gradient: 2W(W^T y - x*)
            x_proj = self.W.T @ y
            grad += 2 * self.W @ (x_proj - x_star)
            
            # L1 subgradient
            grad += self.lambda_1 * np.sign(y)
            
            # Geometric gradient (via chain rule)
            grad += self.lambda_3 * self._compute_R_geo_gradient(y)
            
            return grad
        
        # BFGS optimization with box constraints for feasibility
        result = minimize(
            objective,
            y_init,
            method='L-BFGS-B',
            jac=gradient,
            options={'maxiter': max_iters, 'ftol': 1e-6}
        )
        
        return result.x
    
    def _compute_R_geo(self, y: np.ndarray) -> float:
        """Geometric regularization: E(W^T y) - log(ρ_coset)"""
        x = self.W.T @ y
        
        neighbors, indices = self.lattice.nearest_neighbors(x)
        rho = self.lattice.coset_density(indices)
        
        # Compute energy
        energy = 0.0
        for term in self.energy_terms:
            energy += term.compute(x, neighbors)
        
        # Penalize low coset density
        return energy - np.log(rho + 1e-10)
    
    def _compute_R_geo_gradient(self, y: np.ndarray) -> np.ndarray:
        """Chain rule: ∇_y R_geo = W · ∇_x R_geo"""
        x = self.W.T @ y
        
        neighbors, indices = self.lattice.nearest_neighbors(x)
        
        # Approximate gradient in x-space
        grad_x = np.zeros(8)
        eps = 1e-3
        
        R0 = self._compute_R_geo(y)
        
        for i in range(8):
            y_plus = y.copy()
            y_plus += eps * self.W[:, i]  # Perturb along W basis
            
            R_plus = self._compute_R_geo(y_plus)
            grad_x[i] = (R_plus - R0) / eps
        
        # Map back to y-space
        return self.W @ grad_x
