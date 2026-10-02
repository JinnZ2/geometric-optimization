#!/usr/bin/env python3
"""
budget_audit.py -- what "equal budget" meant in the C20 run, and the per-cell
table with adjusted p.  Re-runs the C20 grid with the same seeds, counting
every call to the landscape function per arm, and records GAS trajectories
so the channel-visit rate can be read against the channel's volume.

Reports, per cell: hits per arm out of 32, Fisher p and Bonferroni-adjusted
p (x12), GAS iterations actually used, landscape-function calls per GAS
iteration and per arm, distinct points evaluated per arm, and the fraction of
runs in which each arm ever PUT A POINT IN THE CHANNEL (theta within 2 s_w of
the wall, psi < w).  Also prints the S^6 cone fraction of the channel per
width, which is the volume rate a uniform sampler lands on it at.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import numpy as np                                            # noqa: E402
from scipy.integrate import quad                              # noqa: E402
from scipy.stats import fisher_exact                          # noqa: E402

import channel_landscape as CL                                # noqa: E402
from gas.solver import GeometricAnnealingSolver, GASParams    # noqa: E402
from validation import experiments as X                       # noqa: E402

WIDTHS, BUDGETS, SEEDS = [0.1, 0.2, 0.4, 0.8, 1.6], [500, 2000], 32
NCELLS = len(WIDTHS) * len(BUDGETS)

_calls = {"n": 0}
_orig = CL.energy_and_grad


def counted(x, w, wall=True):
    _calls["n"] += 1
    return _orig(x, w, wall)


CL.energy_and_grad = counted          # the terms look the name up in the module at call time


def in_channel(x, w):
    n = np.linalg.norm(x); xh = x / n
    theta = np.arccos(np.clip(xh[0], -1, 1)); st = np.sin(theta)
    if st < 1e-9:
        return False
    psi = np.arccos(np.clip(xh[1] / st, -1, 1))
    return abs(theta - CL.THETA_W) <= 2 * CL.S_W and psi < w


def cone_fraction(w):
    """Fraction of S^6 directions within angle w of a fixed direction (azimuthal gap)."""
    num = quad(lambda t: np.sin(t) ** 5, 0, w)[0]
    den = quad(lambda t: np.sin(t) ** 5, 0, np.pi)[0]
    return num / den


def main():
    L = X.lattice()
    rows = []
    print("per-cell table (C20 grid, same seeds), with budget accounting")
    print("  %-6s %-5s | GAS hits  rnd hits  Fisher p   adj p (x12) | GAS iters  calls/iter  GAS calls  rnd calls | GAS pts  rnd pts | in-channel GAS  rnd | cone frac"
          % ("width", "B"))
    for w in WIDTHS:
        cf = cone_fraction(w)
        for B in BUDGETS:
            gh = rh = 0; iters = []; gcalls = []; rcalls = []; gpts = []; g_ch = r_ch = 0
            for seed in range(SEEDS):
                visited = []
                _calls["n"] = 0
                s = GeometricAnnealingSolver(L, CL.terms(w, True), GASParams(max_iters=B),
                                             rng=np.random.default_rng(seed))
                st = s.optimize(callback=lambda state: visited.append(state.x.copy()))
                gcalls.append(_calls["n"]); iters.append(st.iteration)
                gh += CL.angle_to_target(st.best_x) < CL.R_B
                pts = {tuple(np.round(v, 10)) for v in visited}
                gpts.append(len(pts) + 1)                    # + the start point
                g_ch += any(in_channel(v, w) for v in visited)
                _calls["n"] = 0
                P = X.sphere_points(B, seed=1000 + seed)
                E = np.array([CL.energy_and_grad(p, w, True)[0] for p in P])
                rcalls.append(_calls["n"])
                rh += CL.angle_to_target(P[int(np.argmin(E))]) < CL.R_B
                r_ch += any(in_channel(p, w) for p in P)
            p = fisher_exact([[gh, SEEDS - gh], [rh, SEEDS - rh]])[1]
            row = {"width": w, "budget": B, "gas_hits": gh, "random_hits": rh, "seeds": SEEDS,
                   "fisher_p": float(p), "adjusted_p": float(min(1.0, NCELLS * p)),
                   "gas_iters_mean": float(np.mean(iters)), "gas_iters_min": int(min(iters)),
                   "landscape_calls_per_gas_iter": float(np.mean(gcalls) / np.mean(iters)),
                   "gas_landscape_calls_mean": float(np.mean(gcalls)),
                   "random_landscape_calls": int(rcalls[0]),
                   "gas_distinct_points_mean": float(np.mean(gpts)), "random_distinct_points": B,
                   "gas_runs_touching_channel": g_ch, "random_runs_touching_channel": r_ch,
                   "channel_cone_fraction_S6": cf}
            rows.append(row)
            print("  %-6g %-5d | %5d/32  %5d/32  %9.2e  %9.2e  | %8.0f  %9.2f  %9.0f  %9d | %7.0f  %7d | %10d/32  %2d/32 | %.2e"
                  % (w, B, gh, rh, p, row["adjusted_p"], row["gas_iters_mean"], row["landscape_calls_per_gas_iter"],
                     row["gas_landscape_calls_mean"], rcalls[0], row["gas_distinct_points_mean"], B, g_ch, r_ch, cf))
    print("\nbudget accounting")
    print("  one GAS iteration = 1 NEW point (the proposal) + a re-evaluation of the current point + 1 gradient at")
    print("  the current point; with the landscape supplied as two copies that is %.1f landscape-function calls per"
          " iteration, of which 1 distinct new point. Random: 1 call per point, 1 distinct point per call (x2 copies"
          " not used: random calls the function directly)." % rows[0]["landscape_calls_per_gas_iter"])
    print("  So 'equal budget' in C20 meant equal NUMBER OF NEW POINTS (max_iters = B), with GAS also receiving a")
    print("  gradient at each point that random does not get, and GAS spending FEWER new points when its own")
    print("  convergence test stopped it early (iters above).")
    print("\nwidth curve, B=2000: GAS basin-hit rate vs the channel's volume fraction")
    for r in rows:
        if r["budget"] == 2000:
            print("  w=%-4g  cone fraction %.2e   GAS basin hits %2d/32   GAS runs touching channel %2d/32   random runs touching channel %2d/32"
                  % (r["width"], r["channel_cone_fraction_S6"], r["gas_hits"], r["gas_runs_touching_channel"], r["random_runs_touching_channel"]))
    json.dump(rows, open("budget_audit.json", "w"), indent=1)


if __name__ == "__main__":
    main()
