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
def W():
    """(N, 8) projection with orthonormal columns, as the decoder requires."""
    rng = np.random.default_rng(0)
    return np.linalg.qr(rng.standard_normal((16, 8)))[0]


@pytest.fixture
def decoder(lattice, W):
    return ProximalGeometricDecoder(W, lattice, create_energy_suite())


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
        y = np.random.default_rng(1).standard_normal(16)
        val = decoder._compute_R_geo(y)
        assert isinstance(val, (float, np.floating))

    def test_qr_orientation_preserves_n(self):
        """API.py used qr(W.T)[0].T, which collapses (N, 8) to (8, 8)."""
        rng = np.random.default_rng(2)
        raw = rng.standard_normal((20, 8))
        assert np.linalg.qr(raw)[0].shape == (20, 8)
        assert np.linalg.qr(raw.T)[0].T.shape == (8, 8)

    def test_rejects_non_orthonormal_W(self, lattice):
        """The geometric chain rule is only valid when W^T W = I."""
        bad = np.random.default_rng(3).standard_normal((16, 8))
        dec = ProximalGeometricDecoder(bad, lattice, create_energy_suite())
        with pytest.raises(ValueError, match="orthonormal"):
            dec.decode(create_test_point(), max_iters=5)

    def test_geo_gradient_is_nonzero(self, decoder):
        """Revision-1 energy terms made R_geo's gradient ~3e-11, i.e. inert."""
        y = np.random.default_rng(4).standard_normal(16) * 0.1
        assert np.linalg.norm(decoder._compute_R_geo_gradient(y)) > 1e-6
