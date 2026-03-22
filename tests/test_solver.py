"""Tests for the Geometric Annealing Solver."""

import numpy as np
import pytest

from gas.lattice import E8Lattice, create_test_point
from gas.energy_terms import create_energy_suite
from gas.solver import GeometricAnnealingSolver, GASParams, GASState


@pytest.fixture
def lattice():
    return E8Lattice()


@pytest.fixture
def solver(lattice):
    terms = create_energy_suite()
    params = GASParams(max_iters=50)
    return GeometricAnnealingSolver(lattice, terms, params)


class TestGASParams:
    def test_defaults(self):
        p = GASParams()
        assert p.k_neighbors == 24
        assert p.max_iters == 1000
        assert p.eta_0 == 0.1

    def test_custom(self):
        p = GASParams(max_iters=500, rho_min=0.7)
        assert p.max_iters == 500
        assert p.rho_min == 0.7


class TestGASState:
    def test_creation(self):
        x = create_test_point()
        state = GASState(x=x, energy=1.0, rho_coset=0.5, iteration=0)
        assert state.converged is False
        assert state.energy_history == []


class TestSolver:
    def test_phi_rotation_orthogonal(self, solver):
        R = solver.R_phi
        np.testing.assert_allclose(R @ R.T, np.eye(8), atol=1e-12)

    def test_phi_rotation_det_one(self, solver):
        np.testing.assert_allclose(np.linalg.det(solver.R_phi), 1.0, atol=1e-12)

    def test_step_returns_state(self, solver):
        x = create_test_point()
        neighbors, indices = solver.lattice.nearest_neighbors(x)
        rho = solver.lattice.coset_density(indices)
        E = solver._compute_energy(x, neighbors, rho)
        state = GASState(x=x, energy=E, rho_coset=rho, iteration=0)
        new_state = solver.step(state)
        assert isinstance(new_state, GASState)
        assert new_state.iteration == 1

    def test_optimize_runs(self, solver):
        np.random.seed(42)
        result = solver.optimize()
        assert isinstance(result, GASState)
        assert result.iteration > 0
        assert len(result.energy_history) > 0

    def test_optimize_with_init(self, solver):
        x = create_test_point()
        result = solver.optimize(x_init=x)
        assert result.iteration > 0

    def test_optimize_callback(self, solver):
        calls = []
        result = solver.optimize(callback=lambda s: calls.append(s.iteration))
        assert len(calls) == result.iteration

    def test_weights_sum(self, solver):
        weights = solver._compute_weights(0.5)
        # Three weights returned for the three default energy terms
        assert len(weights) == 3
        assert np.all(np.isfinite(weights))
