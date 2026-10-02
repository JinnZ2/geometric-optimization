#!/usr/bin/env python3
"""
experiment_c20.py -- claim C20 through the repo's own claim system.

    python funnel_probe/experiment_c20.py selftest      landscape + gradient checks
    python funnel_probe/experiment_c20.py run           registers the experiment in
                                                        validation.experiments.REGISTRY
                                                        and calls
                                                        `scientific_method run --claim C20`
                                                        (writes claims.json, VALIDATION.md
                                                        through the repo's own tool, and
                                                        c20_results.json into the cwd)

C20 itself is added with `python -m validation.scientific_method add ...`
BEFORE this is run.  Nothing in validation/ or gas/ is edited: the experiment
function enters REGISTRY via the repo's own @experiment decorator at import.

Grid and decision rule are those of PREREGISTRATION.md: widths {0.1, 0.2, 0.4,
0.8, 1.6} rad + NO_WALL control, budgets {500, 2000}, 32 seeds, hit =
angle(best_x, t) < r_b, two-sided Fisher exact per cell, C20 SUPPORTED iff
some cell (control excluded) has p < 0.05/12.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import numpy as np                                   # noqa: E402
from scipy.stats import fisher_exact                 # noqa: E402

from gas.solver import GeometricAnnealingSolver, GASParams   # noqa: E402
from validation import experiments as X              # noqa: E402
from channel_landscape import terms, energy_and_grad, angle_to_target, R_B  # noqa: E402

WIDTHS = [0.1, 0.2, 0.4, 0.8, 1.6]
BUDGETS = [500, 2000]
SEEDS = 32
ALPHA = 0.05 / (len(WIDTHS) * len(BUDGETS))


def gas_run(L, w, wall, budget, seed):
    s = GeometricAnnealingSolver(L, terms(w, wall), GASParams(max_iters=budget),
                                 rng=np.random.default_rng(seed))
    st = s.optimize()                     # x_init=None: uniform start from this rng
    return angle_to_target(st.best_x) < R_B, st.iteration, float(st.best_energy)


def random_run(w, wall, budget, seed):
    pts = X.sphere_points(budget, seed=1000 + seed)   # the repo's own baseline sampler
    E = np.array([energy_and_grad(p, w, wall)[0] for p in pts])
    best = pts[int(np.argmin(E))]
    return angle_to_target(best) < R_B, float(E.min())


def run_grid(verbose=True):
    L = X.lattice()
    cells = []
    for w, wall in [(wd, True) for wd in WIDTHS] + [(0.0, False)]:
        for B in BUDGETS:
            gh = rh = 0; iters = []; ge = []; re = []
            for seed in range(SEEDS):
                h, it, be = gas_run(L, w, wall, B, seed); gh += h; iters.append(it); ge.append(be)
                h2, be2 = random_run(w, wall, B, seed); rh += h2; re.append(be2)
            p = fisher_exact([[gh, SEEDS - gh], [rh, SEEDS - rh]])[1]
            cell = {"width": (w if wall else "NO_WALL"), "budget": B, "gas_hits": gh, "random_hits": rh,
                    "n": SEEDS, "fisher_p": float(p), "gas_iters_mean": float(np.mean(iters)),
                    "gas_iters_min": int(min(iters)), "gas_best_E_mean": float(np.mean(ge)),
                    "random_best_E_mean": float(np.mean(re)),
                    "direction": ("GAS>random" if gh > rh else "GAS<random" if gh < rh else "tie"),
                    "control": not wall}
            cells.append(cell)
            if verbose:
                print("  w=%-8s B=%-5d GAS %2d/%d  random %2d/%d  p=%.2e  %-10s iters %.0f (min %d)" % (
                    cell["width"], B, gh, SEEDS, rh, SEEDS, p, cell["direction"],
                    cell["gas_iters_mean"], cell["gas_iters_min"]))
    return cells


def predictions(cells):
    nw = {c["budget"]: c for c in cells if c["control"]}
    byB = {B: [c for c in cells if not c["control"] and c["budget"] == B] for B in BUDGETS}
    p1 = all(nw[B]["gas_hits"] > nw[B]["random_hits"] for B in BUDGETS)
    p2 = all(all(byB[B][i + 1]["gas_hits"] >= byB[B][i]["gas_hits"] - 2 for i in range(len(WIDTHS) - 1)) for B in BUDGETS)
    # P3: random flat in w -- every width's random count inside the binomial 95% band of the pooled rate
    p3 = True
    for B in BUDGETS:
        rate = np.mean([c["random_hits"] for c in byB[B]]) / SEEDS
        sd = np.sqrt(SEEDS * rate * (1 - rate))
        p3 &= all(abs(c["random_hits"] - SEEDS * rate) <= 2 * sd + 1 for c in byB[B])
    cross = []
    for B in BUDGETS:
        for c in byB[B]:
            if c["gas_hits"] > c["random_hits"]:
                cross.append(c["width"]); break
        else:
            cross.append(None)
    p4 = all(x is None or x >= 0.8 for x in cross)
    return {"P1 NO_WALL: GAS > random at both budgets": p1,
            "P2 GAS hits non-decreasing in w (tolerance 2 of 32)": p2,
            "P3 random hits flat in w (within 2 sd + 1 of the pooled rate)": p3,
            "P4 first width where GAS overtakes random is >= 0.8 rad (or none): %s" % cross: p4}


@X.experiment("thin_channel_vs_random")
def exp_thin_channel_vs_random() -> X.Verdict:
    """C20: on the channel landscape, GAS and uniform random differ in basin hit
    rate at equal budget, in at least one (width, budget) cell at p < 0.05/12."""
    print("\n  C20 grid (PREREGISTRATION.md): widths %s + NO_WALL, budgets %s, %d seeds" % (WIDTHS, BUDGETS, SEEDS))
    cells = run_grid()
    tested = [c for c in cells if not c["control"]]
    sig = [c for c in tested if c["fisher_p"] < ALPHA]
    preds = predictions(cells)
    out = {"cells": cells, "alpha_bonferroni": ALPHA, "significant_cells": sig, "predictions": preds}
    with open("c20_results.json", "w") as fh:
        json.dump(out, fh, indent=1)
    dirs = sorted(set(c["direction"] for c in sig))
    measured = ("%d of %d cells differ at p < %.4f (directions: %s); GAS hits by width at B=2000: %s; "
                "random hits by width at B=2000: %s; NO_WALL control GAS %d vs random %d at B=2000" % (
                    len(sig), len(tested), ALPHA, ",".join(dirs) or "none",
                    [c["gas_hits"] for c in tested if c["budget"] == 2000],
                    [c["random_hits"] for c in tested if c["budget"] == 2000],
                    [c for c in cells if c["control"] and c["budget"] == 2000][0]["gas_hits"],
                    [c for c in cells if c["control"] and c["budget"] == 2000][0]["random_hits"]))
    data = {"n_significant_cells": float(len(sig)), "alpha": ALPHA}
    for c in cells:
        key = "w%s_B%d" % (c["width"], c["budget"])
        data[key + "_gas"] = float(c["gas_hits"]); data[key + "_rnd"] = float(c["random_hits"])
        data[key + "_p"] = float(c["fisher_p"])
    print("  predictions:")
    for k, v in preds.items():
        print("    %s %s" % ("HELD  " if v else "FAILED", k))
    return X.Verdict(len(sig) > 0, measured, data)


def main(argv):
    if not argv or argv[0] == "selftest":
        from channel_landscape import selftest
        return 0 if selftest() else 1
    if argv[0] == "run":
        from validation.scientific_method import main as sm_main
        return sm_main(["run", "--claim", "C20"])
    print(__doc__); return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
