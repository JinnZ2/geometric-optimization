#!/usr/bin/env python3
"""
premise_shelf.py -- why premise_check.py does not converge, read off the landscape.

Written AFTER the full-grid premise check stalled twice (run 1: h=0.02, 78% of
cells unconverged at 10000 steps; run 2: per-cell step, killed at 361x361 after
a 37x37 subgrid showed 73-75% unconverged at the 50000-step cap with the
gradient norm sitting at a constant 3.00e-4).  This file changes no constant and
no rule; it computes where the stalled cells are and why, from `reduced()`.

Two causes, found in this order and recorded in that order.  The first reading
(the shelf below) explained the LOCATION of some stalled cells; the BULK of the
stall turned out to be the theta clip at the antipode holding converged cells
1e-3 short of the pole, where the bowl gradient is exactly A sin(1e-3)/2 =
3.00e-4 > TOL.  That is premise_check's AMENDMENT 2, together with a psi
step bound for the channel's stiffness at small w.  After it, the cells that
stay unconverged are the shelf drifters this file is about.

Outside the channel sector C(u) -> 0 and the reduced energy is

    E(theta) = A (1 + cos theta)/2 + H W(theta)          (W: wall Gaussian at theta_w)

whose theta-derivative has a root INSIDE the wall where the bowl's push toward
the antipode (-A sin theta / 2) balances the wall's push back (H dW/dtheta).
That root is a SHELF: a ring of equilibria in theta, flat in psi except for
the channel's own tail dC/dpsi.  Every start with theta below the wall slides
onto the shelf; from the shelf the only motion is a psi-drift of order the
C-tail, i.e. the basin at (theta_m, psi=0) is reached from the WHOLE shelf
given enough time.  So under gradient flow the basin's attraction set is not
confined to the sector (P3 false in the limit), and the pre-registered
convergence budget cannot reach the limit (drift time >> 1000 time units).

Output: shelf location, the psi-gradient on the shelf per width, the drift
time from psi = pi/2 to the sector edge 3w, and the 37x37 flow census
(unconverged fraction, where they sit) using premise_check's own `descend`.
"""
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np                          # noqa: E402
import channel_landscape_v2 as V2           # noqa: E402
import premise_check as P                   # noqa: E402

NG_DIAG = 37


def dE_dtheta_outside(th):
    E, dEt, dEu = V2.reduced(th, -1.0, 1.0)   # u=-1 (psi=pi): C = exp(-2/w^2) ~ 0 at w=1
    return dEt


def shelf_theta():
    """root of dE/dtheta between the mound and the wall, bisection."""
    lo, hi = V2.THETA_M, V2.THETA_W
    flo = dE_dtheta_outside(lo)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        fm = dE_dtheta_outside(mid)
        if (fm < 0) == (flo < 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main():
    t0 = time.time()
    ts = float(shelf_theta())
    E_s = float(V2.reduced(ts, -1.0, 1.0)[0])
    print("shelf: dE/dtheta = 0 outside the channel at theta* = %.4f (theta_m %.2f, theta_w %.2f); E there %+.4f ; bowl at antipode 0, basin E(t_m) %+.4f"
          % (ts, V2.THETA_M, V2.THETA_W, E_s, float(V2.reduced(V2.THETA_M, 1.0, 0.4)[0])))
    print("  sign of dE/dtheta just inside / just outside theta*: %+.2e / %+.2e  (negative->positive = stable in theta)"
          % (dE_dtheta_outside(ts - 0.02), dE_dtheta_outside(ts + 0.02)))
    out = {"theta_shelf": ts, "E_shelf": E_s, "widths": {}}
    print("\n  w     |dE/dpsi| on the shelf at psi=pi/2   psi=2.5     drift time pi/2 -> 3w (time units)   steps at h=0.02")
    for w in P.WIDTHS:
        g = {}
        for psi in (math.pi / 2, 2.5):
            E, dEt, dEu = V2.reduced(ts, math.cos(psi), w)
            g[psi] = abs(float(dEu) * math.sin(psi)) / math.sin(ts) ** 2     # metric psi-velocity
        # drift time: integrate dpsi / |psi'| from 3w to pi/2 on the shelf
        psis = np.linspace(min(3 * w, math.pi / 2), math.pi / 2, 4001)
        E, dEt, dEu = V2.reduced(np.full_like(psis, ts), np.cos(psis), w)
        vel = np.abs(dEu * np.sin(psis)) / math.sin(ts) ** 2
        tdrift = float(np.trapezoid(1.0 / np.maximum(vel, 1e-300), psis)) if 3 * w < math.pi / 2 else 0.0
        print("  %-4s  %.3e                      %.3e   %.3e                         %.1e" % (w, g[math.pi / 2], g[2.5], tdrift, tdrift / 0.02))
        out["widths"][str(w)] = {"grad_psi_half_pi": g[math.pi / 2], "grad_psi_2p5": g[2.5], "drift_time": tdrift}
    print("  (premise_check budget: 50000 steps = 1000 time units at h=0.02, 200 within the wall band)")

    # flow census on the coarse grid with premise_check's own numerics
    th1 = np.linspace(0.005, math.pi - 0.005, NG_DIAG); ps1 = np.linspace(0, math.pi, NG_DIAG)
    TH, PS = np.meshgrid(th1, ps1, indexing="ij"); TH, PS = TH.ravel(), PS.ravel()
    print("\n  flow census, %d x %d subgrid, premise_check.descend (h %.2f / %.3f, cap %d, tol %.0e):" % (NG_DIAG, NG_DIAG, P.H_STEP, P.H_WALL, P.MAX_STEPS, P.TOL))
    print("  w     unconverged    on shelf |theta-theta*|<0.03   median |grad|   at antipode   near t_m (angle<0.15)")
    for w in P.WIDTHS:
        tt, pp, active, used = P.descend(TH, PS, w)
        g = P.gnorm(tt, pp, w)
        on_shelf = int((active & (np.abs(tt - ts) < 0.03)).sum())
        anti = int((tt > math.pi - 0.05).sum())
        near = int((P.angle_to_tm(tt, pp) < V2.R_H).sum())
        print("  %-4s  %4d/%d      %4d                          %.2e       %4d          %4d"
              % (w, int(active.sum()), TH.size, on_shelf, float(np.median(g[active])) if active.any() else 0.0, anti, near))
        out["widths"][str(w)].update({"cells": int(TH.size), "unconverged": int(active.sum()), "on_shelf": on_shelf,
                                      "antipode": anti, "near_tm": near})
        sys.stdout.flush()
    tot_unconv = sum(v["unconverged"] for v in out["widths"].values())
    tot_shelf = sum(v["on_shelf"] for v in out["widths"].values())
    print("\nREADING (computed): %d unconverged cells over the five widths, %d of them on the shelf ring theta* = %.3f." % (tot_unconv, tot_shelf, ts))
    print("  The shelf is a landscape property (the bowl's outward push and the wall's inward push cancel at theta*).")
    print("  Along psi the shelf is flat to the channel's own tail; the unconverged cells are the ones that tail is")
    print("  still pulling toward the sector edge 3w at a gradient between TOL and 1e-3. In the t -> inf limit they enter")
    print("  the basin, so the basin's gradient-flow attraction set is the sector PLUS whatever part of the shelf")
    print("  the tail reaches; P3 (basin cells inside psi < 3w) is evaluated by premise_check on the CONVERGED cells only,")
    print("  and the unconverged count is printed beside it so the limit is not read as a pass. %.0f s" % (time.time() - t0))
    json.dump(out, open("premise_shelf.json", "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
