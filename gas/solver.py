"""
Geometric Annealing Solver (GAS) - Core Algorithm.

Revision 2. Defects corrected relative to revision 1 (kept in
``legacy/gas/solver.py``), each of which is now pinned by an experiment in
``validation/``:

* ``step()`` scored the proposal ``x_prop`` against the *current* point's
  neighbour set. Combined with revision-1 energy terms that ignored ``x``
  entirely, this made ``E_prop == E_current`` identically, so ``delta_E`` was
  always 0, every move was accepted, and the recorded energy never changed --
  1000 iterations produced a single distinct energy value.
* ``_compute_weights`` returned exactly 3 weights regardless of how many terms
  were supplied, so ``zip()`` silently discarded terms 4-7 and ``--all-terms``
  was a no-op.
* Neither ``sigma_t`` nor ``T_t`` depended on the iteration index, so nothing
  actually annealed; ``sigma_t`` stayed within [0.0452, 0.0500] forever.
* The solver returned its *last* state rather than its *best*, discarding the
  minimum an annealer exists to find.
* Sampling used the global ``np.random`` state, so seeding was process-global.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional, Callable, TYPE_CHECKING

from .energy_terms import is_exceptional

if TYPE_CHECKING:
    from .lattice import E8Lattice
    from .energy_terms import EnergyTerm


@dataclass
class GASParams:
    """Hyperparameters for GAS"""

    k_neighbors: int = 24
    eta_0: float = 0.1  # Initial perturbation
    gamma: float = 2.0  # Annealing rate
    beta: float = 5.0  # Sigmoid steepness
    rho_0: float = 0.5  # Coset density threshold
    # Scaled so that a typical gradient step (|grad| ~ 0.05) is comparable to
    # sigma_0, the noise scale. The revision-1 value of 0.01 multiplied a
    # gradient that was identically zero, so it was never exercised.
    alpha_0: float = 1.0  # Initial learning rate
    sigma_0: float = 0.05  # Initial noise scale
    max_iters: int = 1000
    tau_E: float = 1e-4  # Energy convergence tolerance
    rho_min: float = 0.6  # Minimum coset density
    tau_phi: float = 0.05  # phi-alignment tolerance
    tau_anneal: float = 250.0  # Cooling time constant (iterations)
    T_0: float = 0.05  # Initial Metropolis temperature


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
    best_x: Optional[np.ndarray] = None
    best_energy: float = np.inf


class GeometricAnnealingSolver:
    """Main GAS optimization loop"""

    def __init__(
        self,
        lattice: "E8Lattice",
        energy_terms: List[EnergyTerm],
        params: GASParams,
        rng: Optional[np.random.Generator] = None,
    ):
        if not energy_terms:
            raise ValueError("at least one energy term is required")
        self.lattice = lattice
        self.energy_terms = energy_terms
        self.params = params
        self.rng = rng if rng is not None else np.random.default_rng()
        self.R_phi = self._construct_phi_rotation()

    def _construct_phi_rotation(self) -> np.ndarray:
        """Construct golden rotation matrix in e1-e8 plane"""
        phi = (1 + np.sqrt(5)) / 2
        norm = np.sqrt(2 + phi)

        R = np.eye(8)
        # Golden rotation in the (e1, e8) plane
        c = 1 / norm
        s = phi / norm

        R[0, 0] = c
        R[0, 7] = s
        R[7, 0] = -s
        R[7, 7] = c

        return R

    def _compute_weights(self, rho: float) -> np.ndarray:
        """Dynamic sigmoid weighting based on coset density.

        Returns one weight per energy term -- revision 1 returned exactly three
        regardless of term count, so ``zip()`` silently dropped the rest.
        """
        t = np.tanh(self.params.beta * (rho - self.params.rho_0))
        w_rational = 0.5 * (1 - t)  # favours low rho
        w_exceptional = 0.5 * (1 + t)  # favours high rho

        weights = np.array(
            [
                w_exceptional if is_exceptional(term) else w_rational
                for term in self.energy_terms
            ]
        )
        # Guard against the revision-1 defect recurring: a weight vector
        # shorter than the term list would let zip() drop terms in silence.
        if len(weights) != len(self.energy_terms):  # pragma: no cover
            raise RuntimeError(
                f"weight/term count mismatch: {len(weights)} vs "
                f"{len(self.energy_terms)}"
            )
        return weights

    def _compute_energy(
        self, x: np.ndarray, neighbors: np.ndarray, rho: float
    ) -> float:
        """Calculate total weighted energy"""
        weights = self._compute_weights(rho)

        total = 0.0
        for term, weight in zip(self.energy_terms, weights):
            total += weight * term.compute(x, neighbors)

        return float(total)

    def _compute_gradient(self, x, neighbors, rho):
        """Compute combined gradient from all energy terms"""
        weights = self._compute_weights(rho)
        gradient = np.zeros(8)

        for term, weight in zip(self.energy_terms, weights):
            gradient += weight * term.gradient(x, neighbors)

        return gradient

    def _evaluate(self, x: np.ndarray):
        """Energy and coset density of ``x`` under its *own* neighbourhood."""
        neighbors, indices = self.lattice.nearest_neighbors(
            x, k=self.params.k_neighbors
        )
        rho = self.lattice.coset_density(indices)
        return neighbors, rho

    def step(self, state: GASState) -> GASState:
        """Execute one GAS iteration"""
        # 1. Get neighbourhood of the current point
        neighbors, rho = self._evaluate(state.x)

        # 2. Adaptive annealing schedule. `cooling` is the piece revision 1
        #    was missing: without it nothing depended on the iteration index.
        cooling = np.exp(-state.iteration / self.params.tau_anneal)
        eta_t = self.params.eta_0 * np.exp(-self.params.gamma * rho) * cooling
        alpha_t = self.params.alpha_0 * cooling
        sigma_t = self.params.sigma_0 * cooling
        T_t = self.params.T_0 * np.exp(-self.params.beta * rho) * cooling

        # 3. Compute gradient (analytic)
        gradient = self._compute_gradient(state.x, neighbors, rho)

        # 4. Compose update: gradient + phi-folding + noise
        x_phi = self.R_phi @ state.x

        x_prop = state.x.copy()
        x_prop += -alpha_t * gradient  # Gradient descent
        x_prop += eta_t * (x_phi - state.x)  # phi-folding bias
        x_prop += sigma_t * self.rng.standard_normal(8)  # Annealing noise

        # Normalize to the norm-sqrt(2) sphere
        x_prop = x_prop / np.linalg.norm(x_prop) * np.sqrt(2)

        # 5. Metropolis acceptance. Both energies are evaluated under each
        #    point's *own* neighbourhood, and under a common weighting, so the
        #    comparison is meaningful; revision 1 scored the proposal against
        #    the current point's neighbours and so always found delta_E == 0.
        neighbors_prop, rho_prop = self._evaluate(x_prop)
        E_current = self._compute_energy(state.x, neighbors, rho)
        E_prop = self._compute_energy(x_prop, neighbors_prop, rho)

        delta_E = E_prop - E_current
        accept = delta_E <= 0 or self.rng.random() < np.exp(
            -delta_E / (T_t + 1e-10)
        )

        if accept:
            x_new, E_new, rho_new = x_prop, E_prop, rho_prop
        else:
            x_new, E_new, rho_new = state.x, E_current, rho

        best_x, best_energy = state.best_x, state.best_energy
        if E_new < best_energy:
            best_x, best_energy = x_new.copy(), E_new

        return GASState(
            x=x_new,
            energy=E_new,
            rho_coset=rho_new,
            iteration=state.iteration + 1,
            energy_history=state.energy_history + [E_new],
            rho_history=state.rho_history + [rho_new],
            best_x=best_x,
            best_energy=best_energy,
        )

    def optimize(
        self, x_init: Optional[np.ndarray] = None, callback: Optional[Callable] = None
    ) -> GASState:
        """Run full GAS optimization"""
        if x_init is None:
            x_init = self.rng.standard_normal(8)
            x_init = x_init / np.linalg.norm(x_init) * np.sqrt(2)

        neighbors, rho_init = self._evaluate(x_init)
        E_init = self._compute_energy(x_init, neighbors, rho_init)

        state = GASState(
            x=x_init,
            energy=E_init,
            rho_coset=rho_init,
            iteration=0,
            energy_history=[E_init],
            rho_history=[rho_init],
            best_x=x_init.copy(),
            best_energy=E_init,
        )

        for t in range(self.params.max_iters):
            state = self.step(state)

            if callback:
                callback(state)

            # Check convergence
            if t > 50:
                recent_energies = state.energy_history[-50:]
                energy_stable = (
                    np.std(recent_energies) / (np.mean(recent_energies) + 1e-10)
                    < self.params.tau_E
                )

                coset_sufficient = state.rho_coset > self.params.rho_min

                if energy_stable and coset_sufficient:
                    state.converged = True
                    break

        return state
