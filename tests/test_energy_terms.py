"""Tests for geometric energy term modules.

The revision-1 suite asserted only shapes, types and finiteness. Every one of
those assertions passed while the terms were constant functions with zero
gradients. The tests below are written to fail on each of those defects: see
``TestFalsification``.
"""

import numpy as np
import pytest

from gas.lattice import E8Lattice, create_test_point
from gas.energy_terms import (
    PHI,
    OctahedralEnergy,
    TetrahedralEnergy,
    GoldenEnergy,
    SquareEnergy,
    HexagonalEnergy,
    DodecahedralEnergy,
    IcosahedralEnergy,
    create_energy_suite,
    is_exceptional,
)

ALL_TERMS = [
    OctahedralEnergy,
    TetrahedralEnergy,
    GoldenEnergy,
    SquareEnergy,
    HexagonalEnergy,
    DodecahedralEnergy,
    IcosahedralEnergy,
]


@pytest.fixture(scope="module")
def lattice():
    return E8Lattice()


@pytest.fixture
def neighbors(lattice):
    nbrs, _ = lattice.nearest_neighbors(create_test_point(), k=24)
    return nbrs


@pytest.fixture
def point():
    return create_test_point()


def sphere_points(n, seed=0):
    rng = np.random.default_rng(seed)
    p = rng.standard_normal((n, 8))
    return p / np.linalg.norm(p, axis=1, keepdims=True) * np.sqrt(2)


class TestEnergyTerms:
    @pytest.mark.parametrize("cls", ALL_TERMS)
    def test_compute_returns_scalar(self, cls, point, neighbors):
        assert isinstance(cls().compute(point, neighbors), (float, np.floating))

    @pytest.mark.parametrize("cls", ALL_TERMS)
    def test_compute_non_negative(self, cls, point, neighbors):
        assert cls().compute(point, neighbors) >= 0.0

    @pytest.mark.parametrize("cls", ALL_TERMS)
    def test_gradient_shape_and_finite(self, cls, point, neighbors):
        grad = cls().gradient(point, neighbors)
        assert grad.shape == (8,)
        assert np.all(np.isfinite(grad))

    @pytest.mark.parametrize("cls", ALL_TERMS)
    def test_undefined_at_origin(self, cls, neighbors):
        term = cls()
        if isinstance(term, GoldenEnergy):
            pytest.skip("GoldenEnergy is defined at the origin")
        with pytest.raises(ValueError):
            term.compute(np.zeros(8), neighbors)


class TestFalsification:
    """One test per defect found in the revision-1 audit.

    Each of these fails against ``legacy/gas/energy_terms.py``.
    """

    @pytest.mark.parametrize("cls", ALL_TERMS)
    def test_energy_depends_on_x(self, cls, neighbors):
        """Revision 1: compute() ignored x entirely."""
        term = cls()
        a = term.compute(create_test_point(1), neighbors)
        b = term.compute(create_test_point(2), neighbors)
        assert a != b, f"{cls.__name__} is independent of x"

    @pytest.mark.parametrize("cls", ALL_TERMS)
    def test_gradient_is_nonzero(self, cls, point, neighbors):
        """Revision 1: gradients were exactly 0.0 for every term."""
        grad = cls().gradient(point, neighbors)
        assert np.linalg.norm(grad) > 1e-8, f"{cls.__name__} has a zero gradient"

    @pytest.mark.parametrize("cls", ALL_TERMS)
    def test_gradient_matches_finite_difference(self, cls, lattice):
        """Code verification: d/dx of the implemented energy."""
        term, eps = cls(), 1e-6
        for p in sphere_points(3, seed=11):
            nb, _ = lattice.nearest_neighbors(p, k=24)
            analytic = term.gradient(p, nb)
            numeric = np.zeros(8)
            for i in range(8):
                a, b = p.copy(), p.copy()
                a[i] += eps
                b[i] -= eps
                numeric[i] = (term.compute(a, nb) - term.compute(b, nb)) / (2 * eps)
            rel = np.linalg.norm(analytic - numeric) / max(
                np.linalg.norm(numeric), 1e-12
            )
            assert rel < 1e-5, f"{cls.__name__}: gradient rel err {rel:.2e}"

    @pytest.mark.parametrize("cls", ALL_TERMS)
    def test_not_constant_over_sphere(self, cls, lattice):
        """Revision 1: TetrahedralEnergy == 1/3, GoldenEnergy == |1-phi|."""
        term = cls()
        vals = []
        for p in sphere_points(60, seed=12):
            nb, _ = lattice.nearest_neighbors(p, k=24)
            vals.append(term.compute(p, nb))
        assert np.ptp(vals) > 1e-6, f"{cls.__name__} is constant over the sphere"

    def test_tetrahedral_is_not_one_third(self, lattice):
        """The exact revision-1 value, pinned so the collapse cannot return."""
        term = TetrahedralEnergy()
        for p in sphere_points(20, seed=13):
            nb, _ = lattice.nearest_neighbors(p, k=24)
            assert abs(term.compute(p, nb) - 1 / 3) > 1e-9

    def test_golden_is_not_phi_minus_one(self, lattice):
        """The exact revision-1 value, pinned so the collapse cannot return."""
        term = GoldenEnergy()
        for p in sphere_points(20, seed=14):
            nb, _ = lattice.nearest_neighbors(p, k=24)
            assert abs(term.compute(p, nb) - abs(1 - PHI)) > 1e-9

    def test_terms_are_distinguishable(self, lattice, point):
        """No two terms compute the same value.

        Not a revision-1 regression: those terms did differ numerically. This
        guards the revision-2 design, where several terms share one formula and
        differ only in `target_cosines`, so a copy-paste slip would silently
        make two of them identical.
        """
        nb, _ = lattice.nearest_neighbors(point, k=24)
        values = [cls().compute(point, nb) for cls in ALL_TERMS]
        assert len(set(np.round(values, 9))) == len(ALL_TERMS)


class TestTargetCosines:
    """The documented ideal angles are the ones actually used."""

    def test_tetrahedral_angle(self):
        assert -1 / 3 in TetrahedralEnergy.target_cosines

    def test_hexagonal_angle(self):
        assert 0.5 in HexagonalEnergy.target_cosines

    def test_dodecahedral_uses_inverse_phi_squared(self):
        """Revision 1's docstring claimed 1/sqrt(5) ~ 1/phi^2; they differ."""
        assert np.isclose(max(DodecahedralEnergy.target_cosines), 1 / PHI**2)
        assert not np.isclose(1 / PHI**2, 1 / np.sqrt(5))

    def test_icosahedral_uses_inverse_sqrt_five(self):
        """cos(63.43 deg) is 1/sqrt(5), not 1/phi as revision 1 commented."""
        assert np.isclose(max(IcosahedralEnergy.target_cosines), 1 / np.sqrt(5))
        assert np.isclose(np.cos(np.radians(63.4349)), 1 / np.sqrt(5), atol=1e-5)


class TestCreateEnergySuite:
    def test_default_three_terms(self):
        terms = create_energy_suite()
        assert len(terms) == 3
        assert isinstance(terms[0], OctahedralEnergy)
        assert isinstance(terms[1], TetrahedralEnergy)
        assert isinstance(terms[2], GoldenEnergy)

    def test_include_all_terms(self):
        assert len(create_energy_suite(include_all=True)) == 7

    def test_family_classification(self):
        assert is_exceptional(GoldenEnergy())
        assert is_exceptional(IcosahedralEnergy())
        assert not is_exceptional(OctahedralEnergy())
        assert not is_exceptional(SquareEnergy())
