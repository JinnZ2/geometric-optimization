"""Tests for geometric energy term modules."""

import numpy as np
import pytest

from gas.lattice import E8Lattice, create_test_point
from gas.energy_terms import (
    OctahedralEnergy,
    TetrahedralEnergy,
    GoldenEnergy,
    SquareEnergy,
    HexagonalEnergy,
    DodecahedralEnergy,
    create_energy_suite,
)


@pytest.fixture
def lattice():
    return E8Lattice()


@pytest.fixture
def neighbors(lattice):
    x = create_test_point()
    nbrs, _ = lattice.nearest_neighbors(x, k=24)
    return nbrs


@pytest.fixture
def point():
    return create_test_point()


class TestEnergyTerms:
    @pytest.mark.parametrize(
        "cls",
        [
            OctahedralEnergy,
            TetrahedralEnergy,
            GoldenEnergy,
            SquareEnergy,
            HexagonalEnergy,
            DodecahedralEnergy,
        ],
    )
    def test_compute_returns_scalar(self, cls, point, neighbors):
        term = cls()
        result = term.compute(point, neighbors)
        assert isinstance(result, (float, np.floating))

    @pytest.mark.parametrize(
        "cls",
        [
            OctahedralEnergy,
            TetrahedralEnergy,
            GoldenEnergy,
            SquareEnergy,
            HexagonalEnergy,
            DodecahedralEnergy,
        ],
    )
    def test_compute_non_negative(self, cls, point, neighbors):
        term = cls()
        result = term.compute(point, neighbors)
        assert result >= 0.0

    @pytest.mark.parametrize(
        "cls",
        [OctahedralEnergy, TetrahedralEnergy, GoldenEnergy],
    )
    def test_gradient_shape(self, cls, point, neighbors):
        term = cls()
        grad = term.gradient(point, neighbors)
        assert grad.shape == (8,)

    @pytest.mark.parametrize(
        "cls",
        [OctahedralEnergy, TetrahedralEnergy, GoldenEnergy],
    )
    def test_gradient_finite(self, cls, point, neighbors):
        term = cls()
        grad = term.gradient(point, neighbors)
        assert np.all(np.isfinite(grad))


class TestCreateEnergySuite:
    def test_default_three_terms(self):
        terms = create_energy_suite()
        assert len(terms) == 3
        assert isinstance(terms[0], OctahedralEnergy)
        assert isinstance(terms[1], TetrahedralEnergy)
        assert isinstance(terms[2], GoldenEnergy)

    def test_include_all_six_terms(self):
        terms = create_energy_suite(include_all=True)
        assert len(terms) == 6
