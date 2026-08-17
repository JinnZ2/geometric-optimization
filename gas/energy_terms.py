"""
Geometric Energy Terms (angular-coordination formulation).

Revision 2. The revision-1 formulation (preserved in ``legacy/gas/``) built a
Gram matrix from the *neighbours alone* and never referenced the state ``x``.
That made every term a function of the discrete neighbour set only:

* ``E(x)`` was piecewise constant, taking <=10 distinct values over the whole
  sphere;
* ``TetrahedralEnergy`` was identically ``1/3`` (lambda_min of a rank-<=8 Gram
  matrix of 24 vectors is always 0);
* ``GoldenEnergy`` was identically ``|1 - phi|`` (it ratioed root *norms*, and
  every E8 root has norm sqrt(2));
* the finite-difference gradient perturbed ``x`` while holding ``neighbors``
  fixed, so it returned exactly zero for every term at every point.

Revision 2 measures the angular coordination *between x and each neighbouring
root*, which is what the terms were always described as measuring. Every term
is now a genuine, continuous function of ``x`` with a nonzero analytic
gradient. See ``validation/`` for the experiments that pin this down.
"""

from __future__ import annotations

import numpy as np
from abc import ABC
from typing import Sequence

PHI = (1 + np.sqrt(5)) / 2

#: Default bandwidth of the proximity kernel, in units of E8 root spacing.
#: Neighbours further than a few multiples of this contribute negligibly, so
#: the energy varies smoothly as the k-nearest-neighbour set turns over.
DEFAULT_BANDWIDTH = 0.75


class EnergyTerm(ABC):
    """Base class for geometric energy terms.

    A term scores how far the *angular coordination* of ``x`` relative to its
    neighbouring E8 roots departs from an ideal polytope signature. Concrete
    subclasses supply that signature via :attr:`target_cosines`.
    """

    #: Ideal values of <x_hat, v_hat> for this coordination geometry.
    target_cosines: Sequence[float] = ()

    def __init__(self, bandwidth: float = DEFAULT_BANDWIDTH):
        self.bandwidth = bandwidth

    # -- kernel ---------------------------------------------------------
    def _kernel(self, x: np.ndarray, neighbors: np.ndarray):
        """Proximity weights w_i and their gradients dw_i/dx."""
        delta = x[None, :] - neighbors  # (k, 8)
        sq = np.sum(delta**2, axis=1)
        h2 = self.bandwidth**2
        w = np.exp(-sq / (2 * h2))
        dw = -(w[:, None] * delta) / h2  # (k, 8)
        return w, dw

    def _cosines(self, x: np.ndarray, neighbors: np.ndarray):
        """Cosines c_i = <x_hat, v_hat_i> and their gradients dc_i/dx."""
        nx = np.linalg.norm(x)
        if nx < 1e-12:
            raise ValueError("energy terms are undefined at the origin")
        x_hat = x / nx
        v_hat = neighbors / np.linalg.norm(neighbors, axis=1, keepdims=True)
        c = v_hat @ x_hat  # (k,)
        # d/dx <x/|x|, v_hat> = (v_hat - c * x_hat) / |x|
        dc = (v_hat - c[:, None] * x_hat[None, :]) / nx
        return c, dc

    def _residual(self, c: np.ndarray):
        """Signed distance from each cosine to its nearest target."""
        if not self.target_cosines:
            raise NotImplementedError(
                f"{type(self).__name__} must define target_cosines "
                f"or override compute()/gradient()"
            )
        targets = np.asarray(self.target_cosines, dtype=float)
        diff = c[:, None] - targets[None, :]
        nearest = np.argmin(np.abs(diff), axis=1)
        return diff[np.arange(len(c)), nearest]

    # -- public API -----------------------------------------------------
    def compute(self, x: np.ndarray, neighbors: np.ndarray) -> float:
        """Kernel-weighted mean squared deviation from the ideal coordination."""
        w, _ = self._kernel(x, neighbors)
        c, _ = self._cosines(x, neighbors)
        r = self._residual(c)
        W = np.sum(w)
        if W < 1e-300:
            return 0.0
        return float(np.sum(w * r**2) / W)

    def gradient(self, x: np.ndarray, neighbors: np.ndarray) -> np.ndarray:
        """Analytic gradient of :meth:`compute` with respect to ``x``.

        Verified against central finite differences by
        ``TestFalsification::test_gradient_matches_finite_difference``.
        """
        w, dw = self._kernel(x, neighbors)
        c, dc = self._cosines(x, neighbors)
        r = self._residual(c)
        W = np.sum(w)
        if W < 1e-300:
            return np.zeros(8)
        s = r**2
        E = np.sum(w * s) / W
        # d/dx [ sum(w*s) / sum(w) ]
        #   = [ sum(dw*(s - E)) + sum(2*w*r*dc) ] / W
        grad = dw * (s - E)[:, None] + 2.0 * (w * r)[:, None] * dc
        total: np.ndarray = np.sum(grad, axis=0) / W
        return total


class OctahedralEnergy(EnergyTerm):
    """Octahedral coordination: neighbours on orthogonal axes (cos in {0, +-1}).

    Revision 1 targeted ``lambda_1(G - I) = 0``, i.e. 24 mutually orthogonal
    unit vectors in R^8 -- structurally impossible, so the term never dropped
    below 8.0.
    """

    target_cosines = (0.0, 1.0, -1.0)


class TetrahedralEnergy(EnergyTerm):
    """Tetrahedral coordination: cos(109.47 deg) = -1/3, plus the axial +1."""

    target_cosines = (-1.0 / 3.0, 1.0)


class SquareEnergy(EnergyTerm):
    """Square coordination: strictly 90 deg, cos = 0 (no axial alignment)."""

    target_cosines = (0.0,)


class HexagonalEnergy(EnergyTerm):
    """Hexagonal close packing: cos(60 deg) = 1/2, cos(120 deg) = -1/2."""

    target_cosines = (0.5, -0.5)


class DodecahedralEnergy(EnergyTerm):
    """Dodecahedral coordination: cos = +-1/phi^2 = +-0.38197.

    Revision 1's docstring asserted ``1/sqrt(5) ~ 1/phi^2``; those differ
    (0.44721 vs 0.38197). The value used here is 1/phi^2, matching the name.
    """

    target_cosines = (1.0 / PHI**2, -1.0 / PHI**2)


class IcosahedralEnergy(EnergyTerm):
    """Icosahedral coordination: adjacent-vertex cos = 1/sqrt(5) = 0.44721.

    Revision 1's comment read ``cos(63.43 deg) = 1/phi ~ 0.618``; cos(63.43 deg)
    is in fact 1/sqrt(5) = 0.44721. The correct value is used here, which also
    distinguishes this term from :class:`DodecahedralEnergy`.
    """

    target_cosines = (1.0 / np.sqrt(5.0), -1.0 / np.sqrt(5.0))


class GoldenEnergy(EnergyTerm):
    """phi-ratio spacing between successive neighbour shells around ``x``.

    Revision 1 formed ratios of the neighbour *norms* ``||v_i||``. Every E8
    root has norm sqrt(2), so every ratio was exactly 1 and the term was the
    constant ``|1 - phi| = 0.618034``. Revision 2 uses the distances *from x*
    to each neighbour, which is what "nested scaling" was meant to capture, and
    scores how close successive shell radii come to the golden ratio.
    """

    def _shell_ratios(self, x: np.ndarray, neighbors: np.ndarray):
        d = np.linalg.norm(x[None, :] - neighbors, axis=1)
        order = np.argsort(d)
        d_sorted = np.maximum(d[order], 1e-12)
        r = d_sorted[1:] / d_sorted[:-1]
        return d_sorted, order, r

    def compute(self, x: np.ndarray, neighbors: np.ndarray) -> float:
        if len(neighbors) < 2:
            return 0.0
        _, _, r = self._shell_ratios(x, neighbors)
        return float(np.mean((r - PHI) ** 2))

    def gradient(self, x: np.ndarray, neighbors: np.ndarray) -> np.ndarray:
        if len(neighbors) < 2:
            return np.zeros(8)
        d_sorted, order, r = self._shell_ratios(x, neighbors)
        v_sorted = neighbors[order]
        # dd_i/dx = (x - v_i) / d_i
        dd = (x[None, :] - v_sorted) / d_sorted[:, None]
        da, db = dd[:-1], dd[1:]
        a, b = d_sorted[:-1], d_sorted[1:]
        # r = b/a  =>  dr/dx = db/a - b*da/a^2
        dr = db / a[:, None] - (b / a**2)[:, None] * da
        total: np.ndarray = 2.0 * np.sum((r - PHI)[:, None] * dr, axis=0) / len(r)
        return total


#: Terms whose ideal geometry is rational (integer/rational cosines). Weighted
#: up at low coset density by :meth:`GeometricAnnealingSolver._compute_weights`.
RATIONAL_TERMS = (
    OctahedralEnergy,
    TetrahedralEnergy,
    SquareEnergy,
    HexagonalEnergy,
)

#: Terms whose ideal geometry involves phi. Weighted up at high coset density.
EXCEPTIONAL_TERMS = (
    GoldenEnergy,
    DodecahedralEnergy,
    IcosahedralEnergy,
)


def is_exceptional(term: EnergyTerm) -> bool:
    """True if ``term`` belongs to the phi (exceptional) family."""
    return isinstance(term, EXCEPTIONAL_TERMS)


def create_energy_suite(include_all: bool = False):
    """Create the standard set of geometric energy terms.

    Args:
        include_all: If True, include square, hexagonal, dodecahedral,
            and icosahedral terms in addition to the core three.

    Returns:
        List of EnergyTerm instances.
    """
    terms = [OctahedralEnergy(), TetrahedralEnergy(), GoldenEnergy()]
    if include_all:
        terms.extend(
            [
                SquareEnergy(),
                HexagonalEnergy(),
                DodecahedralEnergy(),
                IcosahedralEnergy(),
            ]
        )
    return terms
