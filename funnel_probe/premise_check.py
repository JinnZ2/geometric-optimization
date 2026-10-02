#!/usr/bin/env python3
"""
premise_check.py -- PREREGISTRATION_C22.md "Premise check": is the lowest
minimum of channel_landscape_v2 reachable ONLY through the channel?
Exhaustive on the reduced (theta, psi) landscape, which is exact for a
landscape depending on x_hat through (theta, psi) alone (asserted by the
module selftest).  Runs BEFORE any GAS or random run; nothing is tuned here.

    grid      theta in linspace(0.005, pi-0.005, 361), psi in linspace(0, pi, 361)
    flow      theta' = -dE/dtheta, psi' = -(1/sin^2 theta) dE/dpsi, RK4 step 0.02,
              until the metric gradient norm < 1e-6 or 10000 steps
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
NG, H_STEP, MAX_STEPS, TOL = 361, 0.02, 10000, 1e-6
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


def descend(th, ps, w):
    th, ps = th.copy(), ps.copy()
    active = np.ones(th.shape, dtype=bool)
    steps_used = np.zeros(th.shape, dtype=int)
    for k in range(MAX_STEPS):
        idx = np.flatnonzero(active)
        if idx.size == 0:
            break
        t, p = th[idx], ps[idx]
        k1 = flow(t, p, w)
        k2 = flow(*wrap(t + 0.5 * H_STEP * k1[0], p + 0.5 * H_STEP * k1[1]), w)
        k3 = flow(*wrap(t + 0.5 * H_STEP * k2[0], p + 0.5 * H_STEP * k2[1]), w)
        k4 = flow(*wrap(t + H_STEP * k3[0], p + H_STEP * k3[1]), w)
        t2 = t + H_STEP / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        p2 = p + H_STEP / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        t2, p2 = wrap(t2, p2)
        th[idx], ps[idx] = t2, p2
        steps_used[idx] = k + 1
        if k % 25 == 0:
            done = gnorm(t2, p2, w) < TOL
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
    print("premise check on the reduced landscape: %d x %d grid, RK4 h=%.2f, tol %.0e, cap %d steps" % (NG, NG, H_STEP, TOL, MAX_STEPS))
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
