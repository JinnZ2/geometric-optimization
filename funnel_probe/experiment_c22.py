#!/usr/bin/env python3
"""
experiment_c22.py -- claim C22 through the repo's own claim system
(PREREGISTRATION_C22.md).  Budget matched on landscape-function CALLS:
the random arm receives exactly the number of calls the GAS run with the
same seed consumed (both term copies counted).

    python funnel_probe/experiment_c22.py run    registers the experiment in
                                                 validation.experiments.REGISTRY and calls
                                                 `scientific_method run --claim C22`
                                                 (writes c22_results.json into the cwd)

Admitted widths are read from premise_check.json (must sit in the cwd or be
passed as C22_PREMISE env var); a width the premise check excluded is not a
cell.  Nothing in validation/ or gas/ is edited.
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

import channel_landscape_v2 as V2                             # noqa: E402
from gas.solver import GeometricAnnealingSolver, GASParams    # noqa: E402
from validation import experiments as X                       # noqa: E402

BUDGETS = [100, 400]          # GAS max_iters; random gets n_calls(GAS) points
SEEDS = 32
_calls = {"n": 0}
_orig = V2.energy_and_grad


def counted(x, w, wall=True):
    _calls["n"] += 1
    return _orig(x, w, wall)


V2.energy_and_grad = counted      # the terms resolve the name in the module at call time


def admitted_widths():
    path = os.environ.get("C22_PREMISE", "premise_check.json")
    d = json.load(open(path))
    return [float(w) for w in d["admitted"]], float(d["f_h"]), d


def run_cell(L, w, M, seed):
    _calls["n"] = 0
    s = GeometricAnnealingSolver(L, V2.terms(w), GASParams(max_iters=M), rng=np.random.default_rng(seed))
    st = s.optimize()
    n_calls = _calls["n"]; it = st.iteration
    gas_hit = V2.angle_to_min(st.best_x) < V2.R_H
    pts = X.sphere_points(n_calls, seed=1000 + seed)
    E = np.array([_orig(p, w)[0] for p in pts])
    rnd_hit = V2.angle_to_min(pts[int(np.argmin(E))]) < V2.R_H
    return gas_hit, rnd_hit, n_calls, it, float(st.best_energy), float(E.min())


@X.experiment("verified_channel_vs_random")
def exp_verified_channel_vs_random() -> X.Verdict:
    """C22: on the exhaustively verified thin-basin landscape, GAS and uniform random
    differ in hit rate (angle(best_x, t_m) < r_h) at equal landscape-call budget, in at
    least one (width, budget) cell at p < 0.05 / n_cells."""
    widths, f_h, prem = admitted_widths()
    cells_n = len(widths) * len(BUDGETS); alpha = 0.05 / cells_n
    print("\n  C22 grid: admitted widths %s, GAS max_iters %s, %d seeds, alpha %.4f ; f_h = %.3e" % (widths, BUDGETS, SEEDS, alpha, f_h))
    L = X.lattice(); cells = []
    for w in widths:
        f_w = prem["widths"][str(w)]["f"]
        for M in BUDGETS:
            gh = rh = 0; calls = []; iters = []; ge = []; re = []
            for seed in range(SEEDS):
                g, r, nc, it, be, bre = run_cell(L, w, M, seed)
                gh += g; rh += r; calls.append(nc); iters.append(it); ge.append(be); re.append(bre)
            p = float(fisher_exact([[gh, SEEDS - gh], [rh, SEEDS - rh]])[1])
            nbar = float(np.mean(calls))
            exp_fh = 1 - (1 - f_h) ** nbar; exp_fw = 1 - (1 - f_w) ** nbar
            cell = {"width": w, "M": M, "gas_hits": gh, "random_hits": rh, "n": SEEDS, "fisher_p": p,
                    "calls_mean": nbar, "calls_min": int(min(calls)), "calls_max": int(max(calls)),
                    "gas_iters_mean": float(np.mean(iters)), "gas_best_E_mean": float(np.mean(ge)),
                    "random_best_E_mean": float(np.mean(re)), "random_expected_from_f_h": exp_fh,
                    "random_expected_from_f_w_prereg_Q1": exp_fw,
                    "direction": "GAS>random" if gh > rh else "GAS<random" if gh < rh else "tie"}
            cells.append(cell)
            print("  w=%-4s M=%-4d GAS %2d/32 random %2d/32 p=%.2e %-10s | calls %.0f (%d-%d) | E_best GAS %.3f rnd %.3f | rnd expected f_h %.4f f(w) %.4f"
                  % (w, M, gh, rh, p, cell["direction"], nbar, min(calls), max(calls), cell["gas_best_E_mean"], cell["random_best_E_mean"], exp_fh, exp_fw))
    sig = [c for c in cells if c["fisher_p"] < alpha]
    sd = lambda n, p: np.sqrt(n * p * (1 - p))
    q1 = all(abs(c["random_hits"] - SEEDS * c["random_expected_from_f_w_prereg_Q1"]) <= 2 * sd(SEEDS, c["random_expected_from_f_w_prereg_Q1"]) + 1 for c in cells)
    q1c = all(abs(c["random_hits"] - SEEDS * c["random_expected_from_f_h"]) <= 2 * sd(SEEDS, c["random_expected_from_f_h"]) + 1 for c in cells)
    q2 = (not any(c["gas_hits"] > c["random_hits"] and c["fisher_p"] < alpha for c in cells)) and \
         all(c["gas_hits"] <= c["random_hits"] or (c["gas_hits"] == 0 and c["random_hits"] == 0) for c in cells if c["width"] <= 0.4)
    q3 = all(c["gas_hits"] == 0 and c["random_hits"] == 0 for c in cells if c["width"] == 0.1)
    preds = {"Q1 random within 2 sd of 1-(1-f(w))^n  [as pre-registered; see premise_check docstring]": q1,
             "Q1' random within 2 sd of 1-(1-f_h)^n  [corrected comparator, stated before the run]": q1c,
             "Q2 GAS does not exceed random at p<alpha anywhere; GAS<=random or both zero at w<=0.4": q2,
             "Q3 both arms 0/32 at w=0.1, both budgets": q3}
    for k, v in preds.items():
        print("  %s %s" % ("HELD  " if v else "FAILED", k))
    json.dump({"cells": cells, "alpha": alpha, "significant": sig, "predictions": preds, "f_h": f_h,
               "admitted_widths": widths}, open("c22_results.json", "w"), indent=1)
    dirs = sorted(set(c["direction"] for c in sig))
    measured = ("%d of %d cells differ at p < %.4f (directions: %s); GAS hits by (w,M): %s; random hits: %s; "
                "calls matched per seed, mean %s; premise verified exhaustively, admitted widths %s" % (
                    len(sig), len(cells), alpha, ",".join(dirs) or "none",
                    [(c["width"], c["M"], c["gas_hits"]) for c in cells], [c["random_hits"] for c in cells],
                    [round(c["calls_mean"]) for c in cells], widths))
    data = {"n_significant_cells": float(len(sig)), "alpha": alpha, "f_h": f_h}
    for c in cells:
        k = "w%s_M%d" % (c["width"], c["M"])
        data[k + "_gas"] = float(c["gas_hits"]); data[k + "_rnd"] = float(c["random_hits"]); data[k + "_p"] = c["fisher_p"]; data[k + "_calls"] = c["calls_mean"]
    return X.Verdict(len(sig) > 0, measured, data)


def main(argv):
    if argv and argv[0] == "run":
        from validation.scientific_method import main as sm_main
        return sm_main(["run", "--claim", "C22"])
    print(__doc__); return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
