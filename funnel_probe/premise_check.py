#!/usr/bin/env python3
"""
premise_check.py -- PREREGISTRATION_C22.md "Premise check": is the lowest
minimum of channel_landscape_v2 reachable ONLY through the channel?
Exhaustive on the reduced (theta, psi) landscape, which is exact for a
landscape depending on x_hat through (theta, psi) alone (asserted by the
module selftest).  Runs BEFORE any GAS or random run; nothing is tuned here.

    grid      theta in linspace(0.005, pi-0.005, 361), psi in linspace(0, pi, 361)
    flow      theta' = -dE/dtheta, psi' = -(1/sin^2 theta) dE/dpsi, RK4 with a
              per-cell step: 0.02 away from the wall, 0.004 within 0.3 rad of
              theta_w (NUMERICAL AMENDMENT, see below), until the metric gradient
              norm < 1e-6 or 50000 steps

NUMERICAL AMENDMENT (before any C22 run; landscape constants untouched): the
first run used h = 0.02 everywhere and 78% of cells never converged -- the
wall term has curvature ~ H/s_w^2 = 200 in theta, so h*lambda = 4 exceeds
RK4's stability bound (~2.8) and cells near theta_w oscillate forever,
scattering into 173 spurious "minima".  The step is reduced only where the
stiffness is; the time horizon (up to 200 time units) is kept.

NUMERICAL AMENDMENT 2 (before any C22 run; landscape constants untouched): run 2
with amendment 1 was stopped after a 37x37 subgrid showed 73-83% of cells still
"unconverged" at the 50000-step cap with the gradient norm at a constant 3.00e-4
-- which is A sin(1e-3)/2, the bowl gradient AT THE THETA CLIP (pi - 1e-3): those
cells had reached the antipode and the clip held them 1e-3 short of the pole, above
TOL, forever.  A cell at either clip now counts as converged (`at_pole`).  The same
subgrid showed a few cells with gradient norm ~20 at w = 0.1: the channel's psi
curvature (D T S + H W) C / (w^2 sin^2 theta) reaches ~900 there, so h = 0.02 is unstable in
psi exactly as h = 0.02 was in theta at the wall.  `step_size` now bounds h by
H_LAMBDA / lambda from a per-cell curvature estimate.  Cells that remain
unconverged after both amendments are reported with their location
(`premise_shelf.py` reads them: a ring of equilibria inside the wall).
    minima    endpoints clustered at 0.02 rad (merged greedily at 0.03)
    volume    cell weight sin^6(theta) sin^5(psi); f(w) = weighted share whose
              descent ends at the global minimum
    PREMISE   P1 one global min within 0.05 rad of t_m, >= 0.2 below every other
              P2 f(w) <= 0.02     P3 basin cells all have psi < 3w
              P4 0 basin cells with theta < theta_w and psi > 3w

Also computed here, before the run and recorded as a correction to the
pre-registration's Q1: the random arm's hit criterion is angle(best_x, t_m)
< r_h, so its comparator is the S^7 cap fraction f_h at r_h, NOT the basin
fraction f(w) the pre-registration wrote.  Both are printed.
"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np                                  # noqa: E402
from scipy.integrate import quad                    # noqa: E402

import channel_landscape_v2 as V2                   # noqa: E402

WIDTHS = [0.1, 0.2, 0.4, 0.8, 1.6]
NG, H_STEP, H_WALL, WALL_BAND, MAX_STEPS, TOL = 361, 0.02, 0.004, 0.3, 50000, 1e-6
H_LAMBDA = 2.0      # RK4 stability target h * lambda <= 2 against a per-cell curvature bound (AMENDMENT 2)
TH_MIN, TH_MAX = 1e-3, np.pi - 1e-3


def flow(th, ps, w):
    """(theta', psi') of the reduced gradient flow, vectorised."""
    E, dEt, dEu = V2.reduced(th, np.cos(ps), w)
    dEp = -dEu * np.sin(ps)                          # dE/dpsi = dE/du * du/dpsi
    s2 = np.maximum(np.sin(th) ** 2, 1e-4)
    return -dEt, -dEp / s2


def gnorm(th, ps, w):
    E, dEt, dEu = V2.reduced(th, np.cos(ps), w)
    dEp = -dEu * np.sin(ps)
    return np.sqrt(dEt ** 2 + dEp ** 2 / np.maximum(np.sin(th) ** 2, 1e-4))


def wrap(th, ps):
    th = np.clip(th, TH_MIN, TH_MAX)
    ps = np.abs(ps); ps = np.pi - np.abs(np.pi - ps)
    return th, ps


def step_size(t, p, w):
    """per-cell RK4 step: the wall band step, then further reduced where the channel is stiff
    (curvature ~ (D T S + H W) C(u) / (w^2 sin^2 theta) in psi), so that h * lambda <= H_LAMBDA."""
    h = np.where(np.abs(t - V2.THETA_W) < WALL_BAND, H_WALL, H_STEP)
    C = np.exp(-(1.0 - np.cos(p)) / (w * w))
    T = np.exp(-(t - V2.THETA_M) ** 2 / (2 * V2.S_M * V2.S_M)) * (1.0 - np.exp(-t * t / (2 * V2.S_0 * V2.S_0)))
    W = np.exp(-(t - V2.THETA_W) ** 2 / (2 * V2.S_W * V2.S_W))
    # the psi-curvature carries the theta envelopes (D T S + H W): without them the 1/sin^2 blows
    # the bound up at the poles where the channel term is physically absent
    lam = (V2.D * T + V2.H * W) * C / (w * w * np.maximum(np.sin(t) ** 2, 1e-4)) + 1.0
    return np.minimum(h, H_LAMBDA / lam)


def at_pole(t):
    """a cell sitting at the theta clip is at a pole of the sphere; the residual gradient there,
    A sin(1e-3)/2 = 3.0e-4, is the clip's not the landscape's (AMENDMENT 2)."""
    return (t >= TH_MAX - 1e-9) | (t <= TH_MIN + 1e-9)


def descend(th, ps, w):
    th, ps = th.copy(), ps.copy()
    active = np.ones(th.shape, dtype=bool)
    steps_used = np.zeros(th.shape, dtype=int)
    for k in range(MAX_STEPS):
        idx = np.flatnonzero(active)
        if idx.size == 0:
            break
        t, p = th[idx], ps[idx]
        h = step_size(t, p, w)
        k1 = flow(t, p, w)
        k2 = flow(*wrap(t + 0.5 * h * k1[0], p + 0.5 * h * k1[1]), w)
        k3 = flow(*wrap(t + 0.5 * h * k2[0], p + 0.5 * h * k2[1]), w)
        k4 = flow(*wrap(t + h * k3[0], p + h * k3[1]), w)
        t2 = t + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        p2 = p + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        t2, p2 = wrap(t2, p2)
        th[idx], ps[idx] = t2, p2
        steps_used[idx] = k + 1
        if k % 25 == 0:
            done = (gnorm(t2, p2, w) < TOL) | at_pole(t2)
            active[idx[done]] = False
    return th, ps, active, steps_used


def angle_to_tm(th, ps):
    c = np.cos(th) * np.cos(V2.THETA_M) + np.sin(th) * np.sin(V2.THETA_M) * np.cos(ps)
    return np.arccos(np.clip(c, -1, 1))


def cluster(th, ps, w):
    key = np.round(th / 0.02).astype(int) * 100000 + np.round(ps / 0.02).astype(int)
    uniq, inv, cnt = np.unique(key, return_inverse=True, return_counts=True)
    cents = np.array([[th[inv == i].mean(), ps[inv == i].mean()] for i in range(len(uniq))])
    # greedy merge at 0.03 in (theta, psi)
    labels = -np.ones(len(uniq), dtype=int); reps = []
    for i in np.argsort(-cnt):
        for j, (ct, cp) in enumerate(reps):
            if np.hypot(cents[i, 0] - ct, cents[i, 1] - cp) < 0.03:
                labels[i] = j; break
        else:
            labels[i] = len(reps); reps.append((cents[i, 0], cents[i, 1]))
    lab = labels[inv]
    mins = []
    for j, (ct, cp) in enumerate(reps):
        E = float(V2.reduced(ct, np.cos(cp), w)[0])
        mins.append({"theta": float(ct), "psi": float(cp), "E": E, "cells": int((lab == j).sum()),
                     "angle_to_tm": float(angle_to_tm(np.array(ct), np.array(cp)))})
    return lab, mins


def cap_fraction_s7(r):
    num = quad(lambda t: np.sin(t) ** 6, 0, r)[0]; den = quad(lambda t: np.sin(t) ** 6, 0, np.pi)[0]
    return num / den


def main():
    t0 = time.time()
    th1 = np.linspace(0.005, np.pi - 0.005, NG); ps1 = np.linspace(0, np.pi, NG)
    TH, PS = np.meshgrid(th1, ps1, indexing="ij"); TH, PS = TH.ravel(), PS.ravel()
    wgt = np.sin(TH) ** 6 * np.sin(PS) ** 5
    f_h = cap_fraction_s7(V2.R_H)
    print("premise check on the reduced landscape: %d x %d grid, RK4 h=%.2f (%.3f within %.1f of the wall, h*lambda <= %.0f in the channel), tol %.0e or at a pole, cap %d steps" % (NG, NG, H_STEP, H_WALL, WALL_BAND, H_LAMBDA, TOL, MAX_STEPS))
    print("r_h = %.2f -> S^7 cap fraction f_h = %.3e (the RANDOM arm's comparator; Q1 as pre-registered used f(w), see docstring)" % (V2.R_H, f_h))
    print("  predicted random hit rate per run: n=600 -> %.4f ; n=2400 -> %.4f" % (1 - (1 - f_h) ** 600, 1 - (1 - f_h) ** 2400))
    results = {"f_h": f_h, "widths": {}}
    admitted = []
    for w in WIDTHS:
        tt, pp, active, used = descend(TH, PS, w)
        lab, mins = cluster(tt, pp, w)
        order = np.argsort([m["E"] for m in mins])
        gmin = mins[order[0]]; second = mins[order[1]] if len(mins) > 1 else None
        gi = int(order[0])
        basin = lab == gi
        f = float(wgt[basin].sum() / wgt.sum())
        p1 = gmin["angle_to_tm"] < 0.05 and (second is None or second["E"] - gmin["E"] >= 0.2)
        p2 = f <= 0.02
        psi_max = float(PS[basin].max()) if basin.any() else 0.0
        p3 = psi_max < 3 * w
        p4_count = int((basin & (TH < V2.THETA_W) & (PS > 3 * w)).sum())
        p4 = p4_count == 0
        holds = p1 and p2 and p3 and p4
        unconv = int(active.sum())
        print("\n  w = %-4s  minima found: %d  (unconverged cells %d, max steps used %d)" % (w, len(mins), unconv, int(used.max())))
        for m in sorted(mins, key=lambda m: m["E"])[:6]:
            print("     E %+.4f at theta %.3f psi %.3f  cells %6d  angle to t_m %.3f" % (m["E"], m["theta"], m["psi"], m["cells"], m["angle_to_tm"]))
        print("     P1 unique global min near t_m, gap >= 0.2 : %s (gap %.3f)" % (p1, (second["E"] - gmin["E"]) if second else float("inf")))
        print("     P2 basin volume fraction f(w) = %.3e <= 0.02 : %s" % (f, p2))
        print("     P3 basin inside sector psi < 3w: max psi %.3f vs %.3f : %s" % (psi_max, 3 * w, p3))
        print("     P4 cap-outside-sector cells in basin: %d : %s" % (p4_count, p4))
        print("     PREMISE %s" % ("HOLDS" if holds else "FAILS"))
        print("     random-arm expectation from f(w) (pre-reg Q1 as written): n=2400 -> %.4f ; from f_h: %.4f" % (1 - (1 - f) ** 2400, 1 - (1 - f_h) ** 2400))
        if holds:
            admitted.append(w)
        results["widths"][str(w)] = {"minima": mins, "f": f, "psi_max_basin": psi_max, "p1": p1, "p2": p2, "p3": p3,
                                     "p4": p4, "p4_count": p4_count, "holds": holds, "unconverged": unconv,
                                     "basin_cells": int(basin.sum())}
    results["admitted"] = admitted
    print("\nADMITTED widths for C22: %s ; EXCLUDED: %s" % (admitted, [w for w in WIDTHS if w not in admitted]))
    print("landscape %s ; %.0f s" % ("IN CLASS" if admitted else "NOT_IN_CLASS", time.time() - t0))
    json.dump(results, open("premise_check.json", "w"), indent=1)
    return 0 if admitted else 1


if __name__ == "__main__":
    sys.exit(main())
