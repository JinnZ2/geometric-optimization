#!/usr/bin/env python3
"""
experiment_c21.py -- claim C21 through the repo's own claim system.

    python funnel_probe/experiment_c21.py selftest   KDE / sampler / hit-criterion checks on
                                                     constructed point sets (no GAS run)
    python funnel_probe/experiment_c21.py run        registers the experiment in
                                                     validation.experiments.REGISTRY and calls
                                                     `scientific_method run --claim C21`
                                                     (writes c21_results.json into the cwd)

Build, grid and decision rule are those of PREREGISTRATION_C21.md, committed
before C21 was added to claims.json and before this file existed:

    density   the C20 GAS arm at B=2000, same 32 seeds per width, every visited
              x recorded; KDE rho(x) = sum_i exp(-|x-x_i|^2 / (2 h^2)), h = 0.3,
              chordal distance on the radius-sqrt(2) sphere
    sampler   rejection from uniform on the sphere, accept with probability
              rho_0 / (rho(x) + rho_0), rho_0 = 0.01 * max rho over a pilot of
              4096 uniform points; rejected candidates cost KDE calls, not
              energy evaluations
    budget    B = 2000 ACCEPTED samples per run, 32 seeds (rng 3000 + seed);
              uniform arm = the C20 random arm (sphere_points, seed 1000 + seed)
    hits      PRIMARY   channel point: |theta - theta_w| <= 2 s_w and psi < w
              SECONDARY basin: angle(best_x, t) < r_b
    decision  C21 SUPPORTED iff some width differs at p < 0.05/5 = 0.01 on the
              PRIMARY criterion (Fisher two-sided on hits out of 32)

Two implementation details the pre-registration left open, fixed here before
the run and marked [CHOICE]: the pilot rng seed (2999) and the candidate
batch size (512, a speed setting with no effect on the accepted sample
distribution).  Nothing in validation/ or gas/ is edited.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import numpy as np                                            # noqa: E402
from scipy.stats import fisher_exact                          # noqa: E402

import channel_landscape as CL                                # noqa: E402
from gas.solver import GeometricAnnealingSolver, GASParams    # noqa: E402
from validation import experiments as X                       # noqa: E402

WIDTHS = [0.1, 0.2, 0.4, 0.8, 1.6]
B = 2000
SEEDS = 32
H = 0.3
RHO0_FRAC = 0.01
PILOT_N = 4096
PILOT_SEED = 2999          # [CHOICE 1] not fixed by the pre-registration
BATCH = 512                # [CHOICE 2] speed only
ALPHA = 0.05 / len(WIDTHS)
RADIUS = np.sqrt(2.0)


def to_sphere(x):
    x = np.asarray(x, dtype=float)
    return x * (RADIUS / np.linalg.norm(x, axis=-1, keepdims=True))


def in_channel(xs, w):
    """PRIMARY criterion, vectorised: |theta - theta_w| <= 2 s_w and psi < w."""
    xh = to_sphere(xs) / RADIUS
    theta = np.arccos(np.clip(xh[..., 0], -1, 1))
    st = np.sin(theta)
    psi = np.arccos(np.clip(np.where(st > 1e-9, xh[..., 1] / np.where(st > 1e-9, st, 1.0), 2.0), -1, 1))
    return (np.abs(theta - CL.THETA_W) <= 2 * CL.S_W) & (psi < w) & (st > 1e-9)


class SphereKDE:
    """rho(x) = sum_i exp(-|x - x_i|^2 / (2 h^2)) for x, x_i on the radius-sqrt2 sphere.
    On the sphere |x - x_i|^2 = 4 - 2 x.x_i, so rho = exp(-2/h^2) * sum_i exp(x.x_i / h^2):
    one matrix product per batch, exact (no cutoff)."""

    def __init__(self, points, h=H):
        self.P = to_sphere(points).astype(np.float64)
        self.h2 = h * h
        self.pref = np.exp(-2.0 / self.h2)

    def __call__(self, xs):
        xs = to_sphere(np.atleast_2d(xs))
        return self.pref * np.exp(xs @ self.P.T / self.h2).sum(axis=1)


def uniform_sphere(n, rng):
    v = rng.normal(size=(n, 8))
    return to_sphere(v)


def inverse_density_sample(kde, rho0, n_accept, rng, batch=BATCH):
    """Rejection sampler. Returns (accepted points, candidates drawn, kde calls)."""
    out = []; got = 0; cand = 0
    while got < n_accept:
        c = uniform_sphere(batch, rng)
        rho = kde(c)
        acc = rng.random(batch) < rho0 / (rho + rho0)
        keep = c[acc]
        out.append(keep); got += len(keep); cand += batch
    pts = np.concatenate(out)[:n_accept]
    return pts, cand


def gas_density_points(L, w, seeds=SEEDS, budget=B):
    pts = []
    for seed in range(seeds):
        s = GeometricAnnealingSolver(L, CL.terms(w, True), GASParams(max_iters=budget),
                                     rng=np.random.default_rng(seed))
        s.optimize(callback=lambda st: pts.append(st.x.copy()))
    return np.array(pts)


def run_width(L, w, verbose=True):
    P = gas_density_points(L, w)
    kde = SphereKDE(P)
    pilot = uniform_sphere(PILOT_N, np.random.default_rng(PILOT_SEED))
    rho_pilot = kde(pilot)
    rho0 = RHO0_FRAC * float(rho_pilot.max())
    inv_ch = inv_bas = uni_ch = uni_bas = 0
    cand_tot = 0; inv_best = []; uni_best = []
    for seed in range(SEEDS):
        rng = np.random.default_rng(3000 + seed)
        xs, cand = inverse_density_sample(kde, rho0, B, rng); cand_tot += cand
        E = np.array([CL.energy_and_grad(p, w, True)[0] for p in xs])          # B energy evaluations
        best = xs[int(np.argmin(E))]
        inv_ch += bool(in_channel(xs, w).any()); inv_bas += bool(CL.angle_to_target(best) < CL.R_B)
        inv_best.append(float(E.min()))
        us = X.sphere_points(B, seed=1000 + seed)                               # the C20 random arm
        Eu = np.array([CL.energy_and_grad(p, w, True)[0] for p in us])
        bu = us[int(np.argmin(Eu))]
        uni_ch += bool(in_channel(us, w).any()); uni_bas += bool(CL.angle_to_target(bu) < CL.R_B)
        uni_best.append(float(Eu.min()))
    p_ch = float(fisher_exact([[inv_ch, SEEDS - inv_ch], [uni_ch, SEEDS - uni_ch]])[1])
    p_bas = float(fisher_exact([[inv_bas, SEEDS - inv_bas], [uni_bas, SEEDS - uni_bas]])[1])
    # where did the accepted mass go? share of accepted samples on the inner side of the wall
    inner_inv = float(np.mean(np.arccos(np.clip(to_sphere(xs)[:, 0] / RADIUS, -1, 1)) < CL.THETA_W))
    inner_uni = float(np.mean(np.arccos(np.clip(to_sphere(us)[:, 0] / RADIUS, -1, 1)) < CL.THETA_W))
    row = {"width": w, "density_points": int(len(P)), "rho0": rho0, "rho_pilot_max": float(rho_pilot.max()),
           "rho_pilot_median": float(np.median(rho_pilot)),
           "candidates_per_accepted": cand_tot / (B * SEEDS),
           "inv_channel_hits": inv_ch, "uni_channel_hits": uni_ch, "p_channel": p_ch,
           "inv_basin_hits": inv_bas, "uni_basin_hits": uni_bas, "p_basin": p_bas,
           "inv_best_E_mean": float(np.mean(inv_best)), "uni_best_E_mean": float(np.mean(uni_best)),
           "inner_cap_share_last_run_inv": inner_inv, "inner_cap_share_last_run_uni": inner_uni,
           "n": SEEDS, "budget_energy_evals_per_run": B}
    if verbose:
        print("  w=%-4s density pts %6d  rho0 %.3g  cand/acc %.2f | channel inv %2d/32 uni %2d/32 p=%.2e |"
              " basin inv %2d/32 uni %2d/32 p=%.2e | inner-cap share inv %.3f uni %.3f" % (
                  w, len(P), rho0, row["candidates_per_accepted"], inv_ch, uni_ch, p_ch,
                  inv_bas, uni_bas, p_bas, inner_inv, inner_uni))
    return row


@X.experiment("inverse_density_vs_random")
def exp_inverse_density_vs_random() -> X.Verdict:
    """C21: sampling inversely to GAS's visit density finds a channel point at a rate
    different from uniform random, in at least one width at p < 0.01 (PRIMARY)."""
    print("\n  C21 grid (PREREGISTRATION_C21.md): widths %s, B=%d accepted samples, %d seeds, h=%.2f" % (
        WIDTHS, B, SEEDS, H))
    L = X.lattice()
    rows = [run_width(L, w) for w in WIDTHS]
    sig = [r for r in rows if r["p_channel"] < ALPHA]
    r1 = any(r["inv_channel_hits"] > r["uni_channel_hits"] and r["p_channel"] < ALPHA for r in rows if r["width"] <= 0.4)
    r2 = len(sig) == 0
    r3 = all(r["inv_basin_hits"] >= r["uni_basin_hits"] for r in rows)
    preds = {"R1 commitment reading: inverse-density > uniform on channel hits at some w <= 0.4 (p<0.01)": r1,
             "R2 bad-luck reading: no width differs at p < 0.01 on channel hits": r2,
             "R3 basin criterion favours inverse density at every width (the stated confound)": r3}
    for k, v in preds.items():
        print("  %s %s" % ("HELD  " if v else "FAILED", k))
    out = {"rows": rows, "alpha": ALPHA, "significant_widths": [r["width"] for r in sig],
           "predictions": preds, "choices": {"PILOT_SEED": PILOT_SEED, "BATCH": BATCH}}
    with open("c21_results.json", "w") as fh:
        json.dump(out, fh, indent=1)
    dirs = sorted(set("inv>uni" if r["inv_channel_hits"] > r["uni_channel_hits"] else "inv<uni" for r in sig))
    measured = ("PRIMARY (channel point) %d of %d widths differ at p < %.3f (directions: %s); channel hits by width "
                "inverse-density %s vs uniform %s; SECONDARY (basin) inverse-density %s vs uniform %s; "
                "density-building cost (32 GAS runs x 2000 per width) excluded by declaration" % (
                    len(sig), len(rows), ALPHA, ",".join(dirs) or "none",
                    [r["inv_channel_hits"] for r in rows], [r["uni_channel_hits"] for r in rows],
                    [r["inv_basin_hits"] for r in rows], [r["uni_basin_hits"] for r in rows]))
    data = {"n_significant_widths": float(len(sig)), "alpha": ALPHA}
    for r in rows:
        k = "w%s" % r["width"]
        for f in ("inv_channel_hits", "uni_channel_hits", "p_channel", "inv_basin_hits", "uni_basin_hits",
                  "p_basin", "rho0", "candidates_per_accepted"):
            data[k + "_" + f] = float(r[f])
    return X.Verdict(len(sig) > 0, measured, data)


def selftest():
    ok = True
    def check(name, cond):
        nonlocal ok; ok &= bool(cond); print(("  PASS " if cond else "  FAIL ") + name)
    print("SELFTEST experiment_c21")
    rng = np.random.default_rng(0)
    P = uniform_sphere(500, rng)
    kde = SphereKDE(P)
    x = P[:3]
    brute = np.array([np.exp(-((xi - P) ** 2).sum(1) / (2 * H * H)).sum() for xi in x])
    check("KDE matrix form == brute-force chordal KDE (rel err < 1e-9)", np.allclose(kde(x), brute, rtol=1e-9))
    far = -P[0]                                           # antipode of a density point
    check("rho at a density point > rho at its antipode", kde(P[0])[0] > kde(far)[0])
    # sampler: with a density concentrated on one cap, accepted mass should leave that cap
    cap = to_sphere(np.c_[np.ones(2000) * 3, rng.normal(size=(2000, 7)) * 0.3])
    k2 = SphereKDE(cap)
    pil = uniform_sphere(PILOT_N, np.random.default_rng(PILOT_SEED))
    rho0 = RHO0_FRAC * k2(pil).max()
    acc, cand = inverse_density_sample(k2, rho0, 1000, np.random.default_rng(1))
    uni = uniform_sphere(1000, np.random.default_rng(2))
    share_acc = np.mean(acc[:, 0] / RADIUS > 0.7); share_uni = np.mean(uni[:, 0] / RADIUS > 0.7)
    check("accepted share in the dense cap (%.3f) < uniform share (%.3f)" % (share_acc, share_uni), share_acc < share_uni)
    check("exactly n_accept points returned", len(acc) == 1000)
    check("sampler can also return a uniform-rate sample when rho0 >> rho (prob -> 1)",
          inverse_density_sample(k2, 1e12, 300, np.random.default_rng(3))[1] == BATCH)
    # channel criterion on constructed points
    th, ps = CL.THETA_W, 0.05
    inside = np.array([np.cos(th), np.sin(th) * np.cos(ps), np.sin(th) * np.sin(ps), 0, 0, 0, 0, 0]) * RADIUS
    outside = np.array([np.cos(th), np.sin(th) * np.cos(1.0), np.sin(th) * np.sin(1.0), 0, 0, 0, 0, 0]) * RADIUS
    offband = np.array([np.cos(th + 3 * CL.S_W), np.sin(th + 3 * CL.S_W) * np.cos(ps), np.sin(th + 3 * CL.S_W) * np.sin(ps), 0, 0, 0, 0, 0]) * RADIUS
    check("channel criterion: point in gap -> True", bool(in_channel(inside, 0.4)))
    check("channel criterion: point behind wall (psi=1.0 > w=0.4) -> False", not bool(in_channel(outside, 0.4)))
    check("channel criterion: point outside the theta band -> False", not bool(in_channel(offband, 0.4)))
    check("Fisher two-sided on identical arms p == 1", fisher_exact([[5, 27], [5, 27]])[1] == 1.0)
    print("selftest:", "OK" if ok else "FAILED")
    return ok


def main(argv):
    if not argv or argv[0] == "selftest":
        return 0 if selftest() else 1
    if argv[0] == "run":
        from validation.scientific_method import main as sm_main
        return sm_main(["run", "--claim", "C21"])
    print(__doc__); return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
