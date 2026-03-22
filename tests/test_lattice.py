"""Tests for the E8 lattice module."""

import numpy as np
import pytest

from gas.lattice import E8Lattice, create_test_point


class TestE8Lattice:
    @pytest.fixture
    def lattice(self):
        return E8Lattice()

    def test_root_count(self, lattice):
        """E8 has exactly 240 roots."""
        assert lattice.all_roots.shape == (240, 8)

    def test_d8_root_count(self, lattice):
        """D8 component has 112 roots."""
        assert lattice.d8_roots.shape == (112, 8)

    def test_coset_root_count(self, lattice):
        """Coset component has 128 roots."""
        assert lattice.coset_roots.shape == (128, 8)

    def test_root_norms(self, lattice):
        """All roots have norm sqrt(2)."""
        norms = np.linalg.norm(lattice.all_roots, axis=1)
        np.testing.assert_allclose(norms, np.sqrt(2), atol=1e-12)

    def test_d8_roots_integer(self, lattice):
        """D8 roots have integer coordinates."""
        assert np.allclose(lattice.d8_roots, np.round(lattice.d8_roots))

    def test_coset_roots_half_integer(self, lattice):
        """Coset roots have half-integer coordinates (±0.5)."""
        np.testing.assert_allclose(np.abs(lattice.coset_roots), 0.5)

    def test_coset_even_parity(self, lattice):
        """Coset roots have even number of negative components."""
        for root in lattice.coset_roots:
            assert np.sum(root < 0) % 2 == 0

    def test_nearest_neighbors_returns_correct_k(self, lattice):
        x = create_test_point()
        neighbors, indices = lattice.nearest_neighbors(x, k=10)
        assert neighbors.shape == (10, 8)
        assert indices.shape == (10,)

    def test_coset_density_range(self, lattice):
        x = create_test_point()
        _, indices = lattice.nearest_neighbors(x, k=24)
        rho = lattice.coset_density(indices)
        assert 0.0 <= rho <= 1.0

    def test_is_coset_mask(self, lattice):
        """First 112 are D8, last 128 are coset."""
        assert not lattice.is_coset[:112].any()
        assert lattice.is_coset[112:].all()


class TestCreateTestPoint:
    def test_returns_8d(self):
        x = create_test_point()
        assert x.shape == (8,)

    def test_norm_sqrt2(self):
        x = create_test_point()
        np.testing.assert_allclose(np.linalg.norm(x), np.sqrt(2), atol=1e-12)

    def test_reproducible(self):
        x1 = create_test_point(seed=123)
        x2 = create_test_point(seed=123)
        np.testing.assert_array_equal(x1, x2)

    def test_different_seeds(self):
        x1 = create_test_point(seed=1)
        x2 = create_test_point(seed=2)
        assert not np.allclose(x1, x2)
