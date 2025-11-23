"""
Geometric Annealing Solver (GAS) - Core Algorithm
"""
import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional, Callable

@dataclass
class GASParams:
    """Hyperparameters for GAS"""
    k_neighbors: int = 24
    eta_0: float = 0.1          # Initial perturbation
    gamma: float = 2.0          # Annealing rate
    beta: float = 5.0           # Sigmoid steepness
    rho_0: float = 0.5          # Coset density threshold
    alpha_0: float = 0.01       # Initial learning rate
    sigma_0: float = 0.05       # Initial noise scale
    max_iters: int = 1000
    tau_E: float = 1e-4        # Energy convergence tolerance
    rho_min: float = 0.6       # Minimum coset density
    tau_phi: float = 0.05      # φ-alignment tolerance


@dataclass
class GASState:
    """Current state of GAS optimization"""
    x: np.ndarray
    energy: float
    rho_coset: float
    iteration: int
    converged: bool = False
    energy_history: List[float] = field(default_factory=list)
    rho_history: List[float] = field(default_factory=list)


class GeometricAnnealingSolver:
    """Main GAS optimization loop"""
    
    def __init__(self, 
                 lattice: 'E8Lattice',
                 energy_terms: List[EnergyTerm],
                 params: GASParams):
        self.lattice = lattice
        self.energy_terms = energy_terms
        self.params = params
        self.R_phi = self._construct_phi_rotation()
    
    def _construct_phi_rotation(self) -> np.ndarray:
        """Construct golden rotation matrix in e₁-e₈ plane"""
        phi = (1 + np.sqrt(5)) / 2
        norm = np.sqrt(2 + phi)
        
        R = np.eye(8)
        # Golden rotation in the (e₁, e₈) plane
        c = 1 / norm
        s = phi / norm
        
        R[0, 0] = c
        R[0, 7] = s
        R[7, 0] = -s
        R[7, 7] = c
        
        return R
    
    def _compute_energy(self, 
                       x: np.ndarray, 
                       neighbors: np.ndarray,
                       rho: float) -> float:
        """Calculate total weighted energy"""
        weights = self._compute_weights(rho)
        
        total = 0.0
        for term, weight in zip(self.energy_terms, weights):
            total += weight * term.compute(x, neighbors)
        
        return total
    
    def _compute_weights(self, rho: float) -> np.ndarray:
        """Dynamic sigmoid weighting based on coset density"""
        # Rational family (oct, tet): favor low rho
        w_rational = 0.5 * (1 - np.tanh(self.params.beta * 
                                       (rho - self.params.rho_0)))
        
        # Exceptional family (phi): favor high rho
        w_exceptional = 0.5 * (1 + np.tanh(self.params.beta * 
                                          (rho - self.params.rho_0)))
        
        # Return weights for [oct, tet, phi, ...] terms
        return np.array([w_rational, w_rational, w_exceptional])
    
    def step(self, state: GASState) -> GASState:
        """Execute one GAS iteration"""
        # 1. Get neighborhood
        neighbors, indices = self.lattice.nearest_neighbors(
            state.x, k=self.params.k_neighbors)
        rho = self.lattice.coset_density(indices)
        
        # 2. Adaptive annealing schedule
        eta_t = self.params.eta_0 * np.exp(-self.params.gamma * rho)
        alpha_t = self.params.alpha_0 / (1 + 0.01 * state.iteration)
        sigma_t = self.params.sigma_0 * np.exp(-eta_t)
        
        # 3. Compute gradient (finite difference)
        gradient = self._compute_gradient(state.x, neighbors, rho)
        
        # 4. Compose update: gradient + φ-folding + noise
        x_phi = self.R_phi @ state.x
        
        x_prop = state.x.copy()
        x_prop += -alpha_t * gradient           # Gradient descent
        x_prop += eta_t * (x_phi - state.x)     # φ-folding bias
        x_prop += sigma_t * np.random.randn(8)  # Annealing noise
        
        # Normalize to S⁷
        x_prop = x_prop / np.linalg.norm(x_prop) * np.sqrt(2)
        
        # 5. Metropolis acceptance
        E_current = state.energy
        E_prop = self._compute_energy(x_prop, neighbors, rho)
        
        T_t = self.params.eta_0 * np.exp(-self.params.beta * rho)
        delta_E = E_prop - E_current
        
        accept = (delta_E <= 0 or 
                 np.random.rand() < np.exp(-delta_E / (T_t + 1e-10)))
        
        if accept:
            new_state = GASState(
                x=x_prop,
                energy=E_prop,
                rho_coset=rho,
                iteration=state.iteration + 1,
                energy_history=state.energy_history + [E_prop],
                rho_history=state.rho_history + [rho]
            )
        else:
            new_state = GASState(
                x=state.x,
                energy=E_current,
                rho_coset=rho,
                iteration=state.iteration + 1,
                energy_history=state.energy_history + [E_current],
                rho_history=state.rho_history + [rho]
            )
        
        return new_state
    
    def _compute_gradient(self, x, neighbors, rho):
        """Compute combined gradient from all energy terms"""
        weights = self._compute_weights(rho)
        gradient = np.zeros(8)
        
        for term, weight in zip(self.energy_terms, weights):
            gradient += weight * term.gradient(x, neighbors)
        
        return gradient
    
    def optimize(self, 
                 x_init: Optional[np.ndarray] = None,
                 callback: Optional[Callable] = None) -> GASState:
        """Run full GAS optimization"""
        if x_init is None:
            x_init = np.random.randn(8)
            x_init = x_init / np.linalg.norm(x_init) * np.sqrt(2)
        
        neighbors, indices = self.lattice.nearest_neighbors(x_init)
        rho_init = self.lattice.coset_density(indices)
        E_init = self._compute_energy(x_init, neighbors, rho_init)
        
        state = GASState(
            x=x_init,
            energy=E_init,
            rho_coset=rho_init,
            iteration=0,
            energy_history=[E_init],
            rho_history=[rho_init]
        )
        
        for t in range(self.params.max_iters):
            state = self.step(state)
            
            if callback:
                callback(state)
            
            # Check convergence
            if t > 50:
                recent_energies = state.energy_history[-50:]
                energy_stable = (np.std(recent_energies) / 
                               (np.mean(recent_energies) + 1e-10) 
                               < self.params.tau_E)
                
                coset_sufficient = state.rho_coset > self.params.rho_min
                
                if energy_stable and coset_sufficient:
                    state.converged = True
                    break
        
        return state
