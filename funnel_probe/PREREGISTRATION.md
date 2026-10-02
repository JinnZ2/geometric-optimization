# PRE-REGISTRATION — funnel_probe/ (geometric-optimization)

WORK ORDER: singular-funnel probes, 4 repos (Kavik via Claude, 2026-10-01).
Source under test: Yanchuk, Wieczorek, Jardon-Kojakhmetov, Alkhayuon,
PRL 137, 147202 (2026); arXiv:2601.02001.
Committed BEFORE any probe code exists and BEFORE claim C20 is registered.
Nothing in `gas/`, `meta_layer/`, `validation/experiments.py` or
`validation/scientific_method.py` is edited. C20 enters `claims.json` only
through `python -m validation.scientific_method add`, and its verdict only
through `... run --claim C20`, with the experiment function registered into
`validation.experiments.REGISTRY` from `funnel_probe/experiment_c20.py` at
import time.

## Outcome enum

    FUNNEL_FOUND         GAS's hit rate depends on channel width (a thin
                         channel narrows the set of starts that reach the
                         basin) and the rate differs from random
    NOT_FOUND_IN_RANGE   no width in the grid where GAS and random differ
                         at the pre-set threshold (grid printed)
    NOT_IN_CLASS         the landscape cannot be made to have a channel that
                         GAS's dynamics see (checked by the no-wall control)
    INCONCLUSIVE         GAS runs hit the iteration cap before the
                         convergence test could decide (reported, not a
                         verdict on C20)

## Claim C20 (registered with the repo's own tool, status UNTESTED, before
## any run)

    "On a landscape whose lowest basin is reached only through a thin
     channel, GAS finds it at a rate different from uniform random at
     equal budget."

    experiment id: thin_channel_vs_random
    C11 stays as recorded. C20's result bears on THIS landscape class only.

## Landscape (ChannelLandscape, an EnergyTerm; depends on direction x_hat only)

    theta = angle(x_hat, t),   t = e_1
    psi   = angle(P_perp x_hat, c),  c = e_2, P_perp = projection off t
            (the azimuth of x around the target axis, in [0, pi])
    E(x)  = 0.25 (1 - cos theta)                              # slope toward t
          + H exp(-(theta - theta_w)^2 / (2 s_w^2)) * G(psi)  # wall
          - D exp(-theta^2 / (2 r_b^2))                        # basin
    G(psi) = 1 - exp(-psi^2 / (2 w^2))                         # channel gap
    D = 1.0, r_b = 0.5 rad, H = 2.0, theta_w = 0.8 rad, s_w = 0.1 rad
    w = channel width (the swept parameter)
    NO_WALL control: H = 0

GAS weights its terms by coset density (`_compute_weights`), rational and
exceptional families summing to 1. The landscape is therefore supplied
TWICE: once as a plain EnergyTerm (rational family) and once as a
GoldenEnergy subclass overriding compute/gradient (exceptional family), so
the weighted total equals E(x) at every rho. Asserted in the selftest.
Analytic gradient checked against central differences (rel err < 1e-5) in
the selftest, as the repo's own C04 does.

## Grid

    widths   w in {0.1, 0.2, 0.4, 0.8, 1.6} rad, plus NO_WALL
    budgets  B in {500, 2000} distinct energy evaluations
             GAS: GASParams(max_iters = B), everything else default;
             random: B uniform points on the norm-sqrt(2) sphere
             (`sphere_points`, the repo's own baseline sampler)
    seeds    32 per cell; GAS rng = default_rng(seed), start x_init = None
             (uniform on the sphere from that rng); random seed = 1000+seed
    hit      angle(best_x, t) < r_b. For random, best_x = argmin E over the
             B points. For GAS, best_x = state.best_x.
    cells    6 widths x 2 budgets = 12

Analytic expectation for random, fixed now: the S^7 cap fraction inside
0.5 rad is ~9.2e-4, so P(hit) ~ 0.37 at B=500 and ~0.84 at B=2000,
independent of w.

## Decision rule (fixed now)

Per cell: two-sided Fisher exact test on (GAS hits, random hits) out of 32
each. C20 SUPPORTED iff at least one cell has p < 0.05/12 (Bonferroni over
the 12 cells); FALSIFIED otherwise. Direction reported per cell. The
NO_WALL control is reported but is NOT one of the 12 cells. GAS iterations
actually used (early convergence) are recorded per run.

Pre-registered predictions:
    P1  NO_WALL: GAS hit rate > random hit rate at both budgets
    P2  GAS hit rate is non-decreasing in w (thin channel traps GAS at
        the wall; the 7-dim azimuthal cone of half-angle w is a fraction
        ~w^6 of directions, so widths <= 0.8 should hold GAS near 0)
    P3  random hit rate is flat in w (within the binomial band)
    P4  the crossing — if any — where GAS overtakes random lies at
        w >= 0.8 rad

## Step 1 — reference + noise (shared across the four repos, identical text)

Reference: `singular_funnel_pitchfork.py` vendored UNCHANGED (sha256 recorded
below). Step 0 result on this machine before any repo work: selftest 4/4 PASS,
STATUS REPRODUCED, mu0=4 eps=0.1 edge log10 x0* = -7.30, slope
d(ln x0*)/d(1/eps) = -1.616. Matches the order's expected values.

Noise script `funnel_noise.py` (stdlib only), Euler-Maruyama, dt = 0.01,
T_max = 600, a = 3, b = 2, eps = 0.1, mu0 = 4.

    model (a)  additive on x, reflecting at x = 0:
               x  <- |x + x(mu - x^2) dt + sigma sqrt(dt) xi|
               mu <- mu + eps(-mu + a x - b) dt
    model (b)  additive on mu, x integrated as y = ln x (deterministic):
               y  <- y + (mu - x^2) dt
               mu <- mu + eps(-mu + a x - b) dt + sigma sqrt(dt) xi

Absorbing exits (declared before any run; differ from the reference's because
the e0 test `x < 1e-8` is not meaningful under x-noise of that size):
    e0  : mu < -0.5 and x < 0.1
    e2  : x > 1.0
    INCONCLUSIVE : neither reached by T_max

Grid:
    starts     mu0 = 4, x0 in {1e-9, 10^-8.5, 1e-8}   (all inside the
               deterministic funnel, whose edge at mu0=4, eps=0.1 is 10^-7.30)
    sigma (a)  10^-11, 10^-10.5, ..., 10^-6      (11 values)
    sigma (b)  10^-4, 10^-3.5, ..., 10^0         (9 values)
    runs       200 per (model, start, sigma); seed = 7919*model_idx
               + 1009*start_idx + 101*sigma_idx + run_idx; random.Random(seed)
    readout    P(e0) with a Wilson 95% interval; sigma_half = the log-
               interpolated sigma at which P(e0) first falls to <= 0.5 of
               P(e0) at the smallest sigma on the grid, per (model, start)

Pre-registered predictions (checked, not tuned):
    P1  at the smallest sigma, P(e0) >= 0.95 for every start (noise
        negligible; reproduces the deterministic funnel)
    P2  P(e0) is non-increasing in sigma up to sampling noise (no rise
        larger than 0.10 between adjacent grid points)
    P3  model (a): sigma_half in [1e-9, 1e-7] (noise competes directly
        with a funnel of width ~5e-8 in x)
    P4  model (b): sigma_half in [1e-2, 1] (noise enters x only through
        the time integral of mu; margin to the edge is ~3.9 nats of ln x)
    P5  sigma_half(b) / sigma_half(a) > 1e4
    P6  within one model, sigma_half across the three starts agrees to
        within a factor of 3 (the width, not the depth, sets survival)

Outcome enum for this step:
    FUNNEL_FOUND         P1 holds and sigma_half lies inside the grid
    NOT_FOUND_IN_RANGE   P(e0) never halves inside the grid (range printed)
    NOT_IN_CLASS         not applicable: this is the paper's own system
    INCONCLUSIVE         > 5% of runs at any grid point used for the
                         decision hit T_max unclassified

vendored reference sha256: 466629a741722183a09221a1594c8f75db7d922ce00eede3afd8e7f01c1f9c46

Scope on every result: deterministic given the seed (GAS is a Metropolis
annealer, so each run is a noisy trajectory with a fixed seed); parameter
set as above; continuous 8-dimensional sphere. A null is a result.
