"""Tests for the proximal-geometric decoder."""

import numpy as np
import pytest

from gas.lattice import E8Lattice, create_test_point
from gas.energy_terms import create_energy_suite
from meta_layer.decoder import ProximalGeometricDecoder


@pytest.fixture
def lattice():
    return E8Lattice()


@pytest.fixture
def decoder(lattice):
    np.random.seed(0)
    W = np.random.randn(16, 8)
    terms = create_energy_suite()
    return ProximalGeometricDecoder(W, lattice, terms)


class TestDecoder:
    def test_decode_returns_correct_dim(self, decoder):
        x_star = create_test_point()
        y = decoder.decode(x_star, max_iters=20)
        assert y.shape == (16,)

    def test_decode_finite(self, decoder):
        x_star = create_test_point()
        y = decoder.decode(x_star, max_iters=20)
        assert np.all(np.isfinite(y))

    def test_geo_reg_returns_scalar(self, decoder):
        y = np.random.randn(16)
        val = decoder._compute_R_geo(y)
        assert isinstance(val, (float, np.floating))
