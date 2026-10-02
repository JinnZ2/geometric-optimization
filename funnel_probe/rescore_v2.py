#!/usr/bin/env python3
"""
rescore_v2.py -- RULES_V2.md applied to EXISTING data (nothing re-run): Rule N (noise band,
bootstrap over seeds) and Rule T (thinness T1, start-selectivity T2) in front of FUNNEL_FOUND.
Usage: python3 funnel_probe/rescore_v2.py   (writes rescore_v2_results.json into the cwd)
"""

import json, os, sys
from math import comb
import numpy as np
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))
DRAWS, SEED, T1_SHARE, ALPHA_T2 = 1000, 0, 0.25, 0.01    # RULES_V2.md; T1_SHARE is [CHOICE 1]


def bootstrap(counts, n, in_F, draws=DRAWS, seed=SEED):
    """Resample each configuration's n seed outcomes: k* ~ Binomial(n, k/n). Returns (|F|*, meanP*) arrays."""
    rng = np.random.default_rng(seed)
    k = np.asarray(counts); p = k / n
    sim = rng.binomial(n, p, size=(draws, len(k))) / n
    return in_F(sim).sum(axis=1), sim.mean(axis=1)


def rule_n(obs, boots, direction):
    """direction: 'nonincreasing' or 'nondecreasing' along the grid order. Returns per-step records."""
    out = []
    for i in range(len(obs) - 1):
        change = obs[i + 1] - obs[i]
        sd = float(np.sqrt(boots[i].var() + boots[i + 1].var()))
        wrong = change > 0 if direction == "nonincreasing" else change < 0
        out.append({"step": i, "change": float(change), "sd": sd, "band": 2 * sd,
                    "violation": bool(wrong and abs(change) > 2 * sd), "wrong_direction": bool(wrong)})
    return out


def two_sided_tails(counts, n, p):
    pmf = np.array([comb(n, j) * p ** j * (1 - p) ** (n - j) for j in range(n + 1)])
    tails = np.array([pmf[pmf <= pmf[j] + 1e-15].sum() for j in range(n + 1)])
    return tails[np.asarray(counts)]


def t_gates(counts, n, in_F):
    k = np.asarray(counts); N = len(k); p = k.mean() / n
    share = float(in_F(k / n).sum() / N)
    t1 = share <= T1_SHARE
    if p <= 0 or p >= 1:
        return {"share": share, "T1": t1, "dispersion": None, "chi2_tail": None, "flagged": 0, "T2": False, "p": p}
    disp = float(k.var() / (n * p * (1 - p)))
    stat = float(((k - n * p) ** 2 / (n * p * (1 - p))).sum()); tail = float(chi2.sf(stat, N - 1))
    flagged = int((two_sided_tails(k, n, p) < 0.05 / N).sum())
    t2 = (disp > 1 and tail < ALPHA_T2) or flagged > 0
    return {"share": share, "T1": t1, "dispersion": disp, "chi2_tail": tail, "flagged": flagged, "T2": t2, "p": p}


def outcome(nonempty, rule_n_ok, persists, gates, t2_evaluable=True):
    if not nonempty:
        return "NOT_FOUND_IN_RANGE (F empty at every drive)"
    if not rule_n_ok:
        return "NOT_FOUND_IN_RANGE (not non-increasing under Rule N)"
    if not persists:
        return "NOT_FOUND_IN_RANGE (does not persist at the slowest drive)"
    if not gates["T1"]:
        return "NOT_FOUND_IN_RANGE, `set not thin` (|F|/N = %.3f > %.2f)" % (gates["share"], T1_SHARE)
    if not t2_evaluable:
        return "FUNNEL_FOUND, qualifier `start-selectivity NOT_EVALUABLE`"
    if not gates["T2"]:
        return "NOT_FOUND_IN_RANGE, `not start-selective` (dispersion %.3f, chi2 tail %.3f, flagged %d)" % (gates["dispersion"], gates["chi2_tail"], gates["flagged"])
    return "FUNNEL_FOUND (V2: thin and start-selective)"


def report_rule_n(name, obs, recs):
    print("  Rule N on %s: %s" % (name, ["%.4g" % o for o in obs]))
    for r in recs:
        print("     step %d: change %+.4g  band +-%.4g  %s" % (r["step"], r["change"], r["band"],
              "VIOLATION" if r["violation"] else ("wrong direction, inside band" if r["wrong_direction"] else "ok")))
    return not any(r["violation"] for r in recs)


def main():
    d = json.load(open(os.path.join(HERE, "samples", "c20_results.json")))
    cells = [c for c in d["cells"] if not c["control"] and c["budget"] == 2000]
    cells.sort(key=lambda c: -c["width"])                      # wide -> thin: the drive slows as w falls
    widths = [c["width"] for c in cells]; hits = [c["gas_hits"] for c in cells]; n = cells[0]["n"]
    print("RULES V2 re-scoring -- GO C20, GAS success set per width at B=2000 (%d seeds = %d starts, one seed each), widths %s" % (n, n, widths))
    # each start is one Bernoulli; bootstrap the 32 outcomes
    rng = np.random.default_rng(SEED)
    boots = [rng.binomial(n, h / n, size=DRAWS) for h in hits]
    okH = report_rule_n("GAS hits (predicted non-increasing as w falls)", hits, rule_n(hits, boots, "nonincreasing"))
    thin = [i for i, h in enumerate(hits) if h > 0]
    i_slow = max(thin) if thin else None
    share = hits[i_slow] / n if thin else None
    g = {"share": share, "T1": share <= T1_SHARE if thin else False, "dispersion": None, "chi2_tail": None, "flagged": 0, "T2": False}
    print("  thinnest width with a non-empty success set: w = %s, share %.3f -> T1 %s ; T2 NOT_EVALUABLE (starts are random seeds, not enumerated)" % (widths[i_slow], share, g["T1"]))
    out_w = outcome(bool(thin), okH, True, g, t2_evaluable=False)
    print("  V2 OUTCOME at w = %s: %s" % (widths[i_slow], out_w))
    for w, h in zip(widths, hits):
        if h == 0:
            print("  V2 OUTCOME at w = %s: NOT_FOUND_IN_RANGE (success set empty, nothing persists)" % w)
    pred_ok = out_w.startswith("FUNNEL_FOUND, qualifier") and widths[i_slow] == 0.4 and all(h == 0 for w, h in zip(widths, hits) if w <= 0.2)
    print("  RULES_V2 prediction for GO (FUNNEL_FOUND at w = 0.4 ONLY with the qualifier; NOT_FOUND at w <= 0.2): %s" % ("HELD" if pred_ok else "FAILED"))
    json.dump({"widths": widths, "hits": hits, "outcome_thinnest": out_w, "prediction_held": bool(pred_ok)},
              open("rescore_v2_results.json", "w"), indent=1)


if __name__ == "__main__":
    main()
