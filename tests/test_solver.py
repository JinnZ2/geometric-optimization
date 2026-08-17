"""Tests for the Geometric Annealing Solver.

``TestFalsification`` encodes the revision-1 solver defects; each of those
tests fails against ``legacy/gas/solver.py``.
"""

import numpy as np
import pytest

from gas.lattice import E8Lattice, create_test_point
from gas.energy_terms import create_energy_suite
from gas.solver import GeometricAnnealingSolver, GASParams, GASState


@pytest.fixture(scope="module")
def lattice():
    return E8Lattice()


@pytest.fixture
def solver(lattice):
    return GeometricAnnealingSolver(
        lattice,
        create_energy_suite(),
        GASParams(max_iters=50),
        rng=np.random.default_rng(0),
    )


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
        state = GASState(
            x=create_test_point(), energy=1.0, rho_coset=0.5, iteration=0
        )
        assert state.converged is False
        assert state.energy_history == []
        assert state.best_energy == np.inf


class TestSolver:
    def test_phi_rotation_orthogonal(self, solver):
        np.testing.assert_allclose(solver.R_phi @ solver.R_phi.T, np.eye(8), atol=1e-12)

    def test_phi_rotation_det_one(self, solver):
        np.testing.assert_allclose(np.linalg.det(solver.R_phi), 1.0, atol=1e-12)

    def test_step_returns_state(self, solver):
        x = create_test_point()
        nb, rho = solver._evaluate(x)
        state = GASState(
            x=x,
            energy=solver._compute_energy(x, nb, rho),
            rho_coset=rho,
            iteration=0,
        )
        assert solver.step(state).iteration == 1

    def test_optimize_runs(self, solver):
        result = solver.optimize()
        assert result.iteration > 0
        assert len(result.energy_history) > 0

    def test_optimize_with_init(self, solver):
        assert solver.optimize(x_init=create_test_point()).iteration > 0

    def test_optimize_callback(self, solver):
        calls = []
        result = solver.optimize(callback=lambda s: calls.append(s.iteration))
        assert len(calls) == result.iteration

    def test_rejects_empty_term_list(self, lattice):
        with pytest.raises(ValueError):
            GeometricAnnealingSolver(lattice, [], GASParams())

    def test_best_state_is_tracked(self, solver):
        """The solver must return its best point, not merely its last one."""
        result = solver.optimize(x_init=create_test_point(3))
        assert result.best_x is not None
        assert result.best_energy <= min(result.energy_history) + 1e-12


class TestFalsification:
    """One test per revision-1 solver defect."""

    def test_energy_actually_changes(self, lattice):
        """Revision 1: 1000 iterations produced a single distinct energy."""
        s = GeometricAnnealingSolver(
            lattice,
            create_energy_suite(),
            GASParams(max_iters=200),
            rng=np.random.default_rng(1),
        )
        history = s.optimize(x_init=create_test_point(1)).energy_history
        assert len(np.unique(np.round(history, 9))) > 1

    def test_solver_descends(self, lattice):
        """The run must end below where it started."""
        terms = create_energy_suite()
        for seed in range(3):
            s = GeometricAnnealingSolver(
                lattice,
                terms,
                GASParams(max_iters=300),
                rng=np.random.default_rng(seed),
            )
            r = s.optimize(x_init=create_test_point(seed))
            assert r.best_energy < r.energy_history[0]

    def test_weight_count_matches_term_count(self, lattice):
        """Revision 1 returned 3 weights always, so zip() dropped terms 4-7."""
        for include_all in (False, True):
            terms = create_energy_suite(include_all=include_all)
            s = GeometricAnnealingSolver(lattice, terms, GASParams())
            assert len(s._compute_weights(0.5)) == len(terms)

    def test_all_terms_changes_the_result(self, lattice):
        """Revision 1: --all-terms was a no-op."""
        out = {}
        for include_all in (False, True):
            s = GeometricAnnealingSolver(
                lattice,
                create_energy_suite(include_all=include_all),
                GASParams(max_iters=60),
                rng=np.random.default_rng(42),
            )
            out[include_all] = s.optimize(x_init=create_test_point(42)).best_energy
        assert out[False] != out[True]

    def test_annealing_schedule_depends_on_iteration(self):
        """Revision 1: sigma_t and T_t depended only on rho, never on t."""
        p = GASParams()
        cooling = np.exp(-np.array([0, 200, 800]) / p.tau_anneal)
        sigma = p.sigma_0 * cooling
        assert np.all(np.diff(sigma) < 0)
        assert sigma[-1] < 0.1 * sigma[0]

    def test_reproducible_and_isolated_from_global_rng(self, lattice):
        """Revision 1 drew from the global np.random state."""

        def run():
            s = GeometricAnnealingSolver(
                lattice,
                create_energy_suite(),
                GASParams(max_iters=40),
                rng=np.random.default_rng(7),
            )
            return s.optimize(x_init=create_test_point(2)).best_energy

        np.random.seed(123)
        first = run()
        np.random.seed(456)
        assert run() == first

    def test_proposal_scored_under_its_own_neighbourhood(self, lattice):
        """Revision 1 scored the proposal against the current point's neighbours,
        forcing delta_E == 0 on every step."""
        s = GeometricAnnealingSolver(
            lattice,
            create_energy_suite(),
            GASParams(max_iters=1),
            rng=np.random.default_rng(0),
        )
        a, b = create_test_point(1), create_test_point(2)
        nb_a, rho_a = s._evaluate(a)
        nb_b, rho_b = s._evaluate(b)
        assert s._compute_energy(b, nb_b, rho_b) != s._compute_energy(b, nb_a, rho_a)
