"""
channel_landscape_v2.py -- a landscape whose lowest minimum's ENTIRE basin of
attraction is a thin channel (PREREGISTRATION_C22.md).  Depends on x_hat only.

    theta = angle(x_hat, e_1)         ct = cos theta
    psi   = azimuth about e_1 from e_2 u  = cos psi = x_hat_2 / sin theta
    E = A (1 + ct)/2                              outer bowl, slope AWAY from the pole
      + H W(theta) (1 - C(u))                     wall ridge with a gap in the sector
      - D C(u) T(theta) S(theta)                  trough inside the sector, deepest at theta_m
    W = exp(-(theta-theta_w)^2/(2 s_w^2))   C = exp(-(1-u)/w^2)
    T = exp(-(theta-theta_m)^2/(2 s_m^2))   S = 1 - exp(-theta^2/(2 s_0^2))

AMENDMENT, made before the premise check and before any run, recorded here:
the pre-registration's formula had no S(theta).  Without it the trough term
is psi-dependent at theta = 0, where psi is undefined, so E has a cone point
at the pole and the gradient is discontinuous there.  S(theta) with
s_0 = 0.1 sends the sector term to zero at the pole and changes nothing at
theta >= 0.4 (S(0.4) = 0.9997).  Constants A, H, theta_w, s_w, D, theta_m,
s_m are exactly the pre-registered ones.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

from gas.energy_terms import EnergyTerm, GoldenEnergy

A, H, THETA_W, S_W, D, THETA_M, S_M, S_0 = 0.6, 2.0, 0.8, 0.1, 1.5, 0.4, 0.35, 0.1
R_H = 0.15
E1 = np.eye(8)[0]
E2 = np.eye(8)[1]
T_M = np.sqrt(2.0) * (np.cos(THETA_M) * E1 + np.sin(THETA_M) * E2)   # the intended global minimum


def reduced(theta, u, w):
    """E(theta, u) and partials dE/dtheta, dE/du.  Vectorised over numpy arrays."""
    theta = np.asarray(theta, dtype=float); u = np.asarray(u, dtype=float)
    ct, st = np.cos(theta), np.sin(theta)
    W = np.exp(-(theta - THETA_W) ** 2 / (2 * S_W * S_W)); dW = -(theta - THETA_W) / (S_W * S_W) * W
    T = np.exp(-(theta - THETA_M) ** 2 / (2 * S_M * S_M)); dT = -(theta - THETA_M) / (S_M * S_M) * T
    S = 1.0 - np.exp(-theta * theta / (2 * S_0 * S_0)); dS = theta / (S_0 * S_0) * np.exp(-theta * theta / (2 * S_0 * S_0))
    C = np.exp(-(1.0 - u) / (w * w)); dC = C / (w * w)
    TS, dTS = T * S, dT * S + T * dS
    E = A * (1.0 + ct) / 2 + H * W * (1.0 - C) - D * C * TS
    dE_dtheta = -A * st / 2 + H * dW * (1.0 - C) - D * C * dTS
    dE_du = -H * W * dC - D * dC * TS
    return E, dE_dtheta, dE_du


def angles(x):
    n = np.linalg.norm(x); xh = x / n
    ct = float(np.clip(xh[0], -1.0, 1.0)); theta = np.arccos(ct)
    st = np.sqrt(max(0.0, 1.0 - ct * ct))
    u = float(np.clip(xh[1] / st, -1.0, 1.0)) if st > 1e-12 else 1.0
    return n, xh, ct, st, theta, u


def energy_and_grad(x, w, wall=True):
    """wall is accepted for interface parity with v1 and ignored (no NO_WALL control in C22)."""
    n, xh, ct, st, theta, u = angles(x)
    E, dE_dtheta, dE_du = reduced(theta, u, w)
    dct = (E1 - ct * xh) / n
    if st > 1e-9:
        dtheta = -dct / st
        dxh1 = (E2 - xh[1] * xh) / n
        du = dxh1 / st + xh[1] * ct / (st ** 3) * dct
        g = float(dE_dtheta) * dtheta + float(dE_du) * du
    else:                                   # at a pole: S -> 0 kills the u-dependence; bowl gradient only
        g = -A / 2 * dct
    return float(E), g


class ChannelV2(EnergyTerm):
    target_cosines = (0.0,)

    def __init__(self, w, wall=True):
        super().__init__(); self.w = float(w)

    def compute(self, x, neighbors):
        return energy_and_grad(x, self.w)[0]

    def gradient(self, x, neighbors):
        return energy_and_grad(x, self.w)[1]


class ChannelV2X(GoldenEnergy):
    def __init__(self, w, wall=True):
        super().__init__(); self.w = float(w)

    def compute(self, x, neighbors):
        return energy_and_grad(x, self.w)[0]

    def gradient(self, x, neighbors):
        return energy_and_grad(x, self.w)[1]


def terms(w):
    return [ChannelV2(w), ChannelV2X(w)]


def angle_to_min(x):
    return float(np.arccos(np.clip(x @ T_M / (np.linalg.norm(x) * np.linalg.norm(T_M)), -1.0, 1.0)))


def selftest():
    from gas.lattice import E8Lattice
    from gas.solver import GeometricAnnealingSolver, GASParams
    ok = True
    def check(name, cond):
        nonlocal ok; ok &= bool(cond); print(("  PASS " if cond else "  FAIL ") + name)
    rng = np.random.default_rng(3); worst = 0.0; worst_red = 0.0
    for _ in range(200):
        x = rng.standard_normal(8); x = x / np.linalg.norm(x) * np.sqrt(2)
        if abs(x[0] / np.sqrt(2)) > 0.995:
            continue
        for w in (0.1, 0.4, 1.6):
            E, g = energy_and_grad(x, w)
            _, _, _, _, th, u = angles(x)
            worst_red = max(worst_red, abs(E - float(reduced(th, u, w)[0])))
            fd = np.zeros(8); h = 1e-6
            for i in range(8):
                xp, xm = x.copy(), x.copy(); xp[i] += h; xm[i] -= h
                fd[i] = (energy_and_grad(xp, w)[0] - energy_and_grad(xm, w)[0]) / (2 * h)
            worst = max(worst, np.linalg.norm(g - fd) / (np.linalg.norm(fd) + 1e-12))
    check("8D energy == reduced E(theta, psi), max |diff| %.1e < 1e-9" % worst_red, worst_red < 1e-9)
    check("analytic gradient vs central FD, rel err %.1e < 1e-5" % worst, worst < 1e-5)
    Em = energy_and_grad(T_M, 0.4)[0]; Ep = energy_and_grad(np.sqrt(2) * E1, 0.4)[0]; Ea = energy_and_grad(-np.sqrt(2) * E1, 0.4)[0]
    check("E(t_m) = %.3f < E(antipode) = %.3f < E(pole) = %.3f" % (Em, Ea, Ep), Em < Ea < Ep)
    off = np.sqrt(2) * (np.cos(THETA_M) * E1 - np.sin(THETA_M) * E2)      # same theta, psi = pi
    check("same theta outside the sector is a hill, E = %.3f > 0" % energy_and_grad(off, 0.4)[0], energy_and_grad(off, 0.4)[0] > 0)
    gp = np.linalg.norm(energy_and_grad(np.sqrt(2) * E1 + 1e-7 * E2, 0.4)[1] - energy_and_grad(np.sqrt(2) * E1 - 1e-7 * E2, 0.4)[1])
    check("gradient continuous across the pole, jump %.1e < 1e-4" % gp, gp < 1e-4)
    L = E8Lattice(); s = GeometricAnnealingSolver(L, terms(0.3), GASParams(), rng=np.random.default_rng(0)); wd = 0.0
    for _ in range(20):
        x = rng.standard_normal(8); x = x / np.linalg.norm(x) * np.sqrt(2)
        for rho in (0.0, 0.3, 0.7, 1.0):
            wd = max(wd, abs(s._compute_energy(x, np.zeros((0, 8)), rho) - energy_and_grad(x, 0.3)[0]))
    check("GAS weighted total == E at rho in {0,.3,.7,1}: max |diff| %.1e" % wd, wd < 1e-9)
    print("selftest:", "OK" if ok else "FAILED"); return ok


if __name__ == "__main__":
    import sys; sys.exit(0 if selftest() else 1)
