"""
channel_landscape.py -- a landscape whose lowest basin is reached only through
a thin channel, as an EnergyTerm for the repo's GAS solver.  Depends on the
direction x_hat only (GAS lives on the norm-sqrt(2) sphere).

    theta = angle(x_hat, t),  t = e_1          (polar angle from the target)
    u     = cos(psi), psi = angle(P_perp x_hat, c),  c = e_2
            (azimuth of x about the target axis; u in [-1, 1])
    E     = 0.25 (1 - cos theta)                             slope toward t
          + H * W(theta) * G(u)                              wall with a gap
          - D * exp(-theta^2 / (2 r_b^2))                    basin
    W     = exp(-(theta - theta_w)^2 / (2 s_w^2))
    G     = 1 - exp(-(1 - u) / w^2)      (~ 1 - exp(-psi^2/(2 w^2)) for small psi;
                                          smooth in u everywhere, unlike a
                                          gaussian in psi, which has a cone
                                          point on the antipodal azimuth)

Constants fixed in PREREGISTRATION.md: D=1, r_b=0.5, H=2, theta_w=0.8, s_w=0.1.
`w` is the swept channel width.  `wall=False` is the NO_WALL control (H=0).

GAS weights its terms by coset density, with the rational and exceptional
families summing to 1, so the landscape is supplied twice (ChannelLandscape is
rational; ChannelLandscapeX subclasses GoldenEnergy and is exceptional) and the
weighted total is E(x) at every rho.  Asserted in selftest().
"""
import numpy as np

from gas.energy_terms import EnergyTerm, GoldenEnergy

D, R_B, H, THETA_W, S_W = 1.0, 0.5, 2.0, 0.8, 0.1
E1 = np.eye(8)[0]
E2 = np.eye(8)[1]


def energy_and_grad(x, w, wall=True):
    n = np.linalg.norm(x)
    xh = x / n
    ct = float(np.clip(xh[0], -1.0, 1.0))
    theta = np.arccos(ct)
    st = np.sqrt(max(0.0, 1.0 - ct * ct))
    dct = (E1 - ct * xh) / n                       # d cos(theta) / dx
    # background slope toward t
    E = 0.25 * (1.0 - ct)
    g = -0.25 * dct
    # basin: d(-D Bs)/dx = -D (theta / (r_b^2 sin theta)) Bs dct, theta/sin -> 1
    Bs = np.exp(-theta * theta / (2 * R_B * R_B))
    ratio = theta / st if st > 1e-12 else 1.0
    E -= D * Bs
    g -= D * (ratio / (R_B * R_B)) * Bs * dct
    if wall:
        W = np.exp(-(theta - THETA_W) ** 2 / (2 * S_W * S_W))
        if st > 1e-9:
            u = float(np.clip(xh[1] / st, -1.0, 1.0))
            G = 1.0 - np.exp(-(1.0 - u) / (w * w))
            dG_du = -np.exp(-(1.0 - u) / (w * w)) / (w * w)   # dG/du < 0: the gap closes as u falls
            dxh1 = (E2 - xh[1] * xh) / n
            du = dxh1 / st + xh[1] * ct / (st ** 3) * dct
            dW = ((theta - THETA_W) / (S_W * S_W * st)) * W * dct
            E += H * W * G
            g += H * (dW * G + W * dG_du * du)
        else:                                      # on the axis: W ~ e^-32, G undefined
            E += H * W
    return float(E), g


class ChannelLandscape(EnergyTerm):
    """Rational-family copy of the channel landscape (ignores neighbours)."""

    target_cosines = (0.0,)   # unused; satisfies the base class contract

    def __init__(self, w, wall=True):
        super().__init__()
        self.w, self.wall = float(w), bool(wall)

    def compute(self, x, neighbors):
        return energy_and_grad(x, self.w, self.wall)[0]

    def gradient(self, x, neighbors):
        return energy_and_grad(x, self.w, self.wall)[1]


class ChannelLandscapeX(GoldenEnergy):
    """Exceptional-family copy (is_exceptional() is an isinstance test)."""

    def __init__(self, w, wall=True):
        super().__init__()
        self.w, self.wall = float(w), bool(wall)

    def compute(self, x, neighbors):
        return energy_and_grad(x, self.w, self.wall)[0]

    def gradient(self, x, neighbors):
        return energy_and_grad(x, self.w, self.wall)[1]


def terms(w, wall=True):
    return [ChannelLandscape(w, wall), ChannelLandscapeX(w, wall)]


def angle_to_target(x):
    return float(np.arccos(np.clip(x[0] / np.linalg.norm(x), -1.0, 1.0)))


def selftest():
    from gas.lattice import E8Lattice
    from gas.solver import GeometricAnnealingSolver, GASParams
    ok = True
    def check(name, cond):
        nonlocal ok
        print(("  PASS " if cond else "  FAIL ") + name); ok &= cond
    rng = np.random.default_rng(3)
    # 1. analytic gradient vs central differences, away from the poles
    worst = 0.0
    for _ in range(40):
        x = rng.standard_normal(8); x = x / np.linalg.norm(x) * np.sqrt(2)
        if abs(x[0] / np.sqrt(2)) > 0.98:
            continue
        for w in (0.1, 0.4, 1.6):
            _, g = energy_and_grad(x, w)
            fd = np.zeros(8); h = 1e-6
            for i in range(8):
                xp, xm = x.copy(), x.copy(); xp[i] += h; xm[i] -= h
                fd[i] = (energy_and_grad(xp, w)[0] - energy_and_grad(xm, w)[0]) / (2 * h)
            worst = max(worst, np.linalg.norm(g - fd) / (np.linalg.norm(fd) + 1e-12))
    check("analytic gradient vs central FD, rel err %.1e < 1e-5" % worst, worst < 1e-5)
    # 2. landscape shape
    t = np.sqrt(2) * E1
    check("global minimum at the target, E(t) = -1", abs(energy_and_grad(t, 0.2)[0] + 1.0) < 1e-12)
    on_wall_gap = np.sqrt(2) * (np.cos(THETA_W) * E1 + np.sin(THETA_W) * E2)       # psi = 0
    on_wall_far = np.sqrt(2) * (np.cos(THETA_W) * E1 - np.sin(THETA_W) * E2)       # psi = pi
    Eg, Ef = energy_and_grad(on_wall_gap, 0.2)[0], energy_and_grad(on_wall_far, 0.2)[0]
    check("wall absent in the channel (psi=0) and present opposite it: %.3f vs %.3f" % (Eg, Ef),
          Eg < 0.1 and Ef > 1.5)
    check("NO_WALL control removes the wall", abs(energy_and_grad(on_wall_far, 0.2, wall=False)[0] - Eg) < 1e-9)
    # 3. GAS weighting: rational + exceptional copies sum to E at any rho
    L = E8Lattice()
    s = GeometricAnnealingSolver(L, terms(0.3), GASParams(), rng=np.random.default_rng(0))
    worst_w = 0.0
    for _ in range(20):
        x = rng.standard_normal(8); x = x / np.linalg.norm(x) * np.sqrt(2)
        nb, idx = L.nearest_neighbors(x, k=24)
        for rho in (0.0, 0.25, 0.5, 0.75, 1.0):
            worst_w = max(worst_w, abs(s._compute_energy(x, nb, rho) - energy_and_grad(x, 0.3)[0]))
    check("solver weighted total equals E(x) at rho in {0,...,1}, max |diff| %.1e" % worst_w, worst_w < 1e-12)
    from gas.energy_terms import is_exceptional
    check("one copy rational, one exceptional", (not is_exceptional(terms(0.3)[0])) and is_exceptional(terms(0.3)[1]))
    return ok


if __name__ == "__main__":
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    sys.exit(0 if selftest() else 1)
