# PRE-REGISTRATION — claim C22 (geometric-optimization)

Order (Kavik via Claude, 2026-10-02): "rebuild the landscape so the basin
is reachable ONLY through the channel; verify that premise by exhaustive
check BEFORE any run (new claim C22)." Committed before
`channel_landscape_v2.py` exists and before C22 is added to claims.json.

## Why C20's landscape did not have the premise

In the C20 landscape the lowest basin is the whole cap theta < theta_w,
and a uniform sampler lands in it at the cap's volume rate from the FAR
side of the wall; the channel was the route for a descending dynamic,
not for a point sampler (budget audit, width curve). So "reachable only
through the channel" was true of gradient flow outside the wall and false
of sampling. The rebuilt landscape must make the lowest minimum's ENTIRE
basin of attraction thin.

## Landscape v2 (declared now; depends on x_hat only)

    theta = angle(x_hat, e_1),  psi = azimuth of x_hat about e_1 from e_2  (as C20)
    ct = cos theta, u = cos psi
    E  = A (1 + ct) / 2                                  outer bowl: 0 at the antipode,
                                                         A at the pole, slope AWAY from the pole
       + H W(theta) (1 - C(u))                           wall ridge with a gap in the channel sector
       - D C(u) T(theta)                                 a trough inside the sector, deepest at theta_m
    W(theta) = exp(-(theta - theta_w)^2 / (2 s_w^2))
    C(u)     = exp(-(1 - u) / w^2)                       channel sector profile (~ exp(-psi^2/(2w^2)))
    T(theta) = exp(-(theta - theta_m)^2 / (2 s_m^2))
    A = 0.6, H = 2.0, theta_w = 0.8, s_w = 0.1, D = 1.5, theta_m = 0.4, s_m = 0.35
    w = channel width (swept): {0.1, 0.2, 0.4, 0.8, 1.6}

Intended structure (to be VERIFIED, not assumed): the global minimum sits
at t_m = (cos theta_m, sin theta_m, 0, ...) inside the cap and inside the
sector; everywhere outside the sector the cap is a HILL (the bowl slopes
away from the pole) so a point landing in the cap outside the sector does
not descend to the minimum; the trough is the only descending route in.

## Premise check (exhaustive on the reduced 2D landscape; BEFORE any run)

The energy depends on (theta, psi) only, so gradient flow on S^7 reduces
exactly to flow on [0, pi] x [0, pi] with metric ds^2 = dtheta^2 +
sin^2(theta) dpsi^2: theta' = -dE/dtheta, psi' = -(1/sin^2 theta) dE/dpsi.

    grid       theta in linspace(0.005, pi - 0.005, 361), psi in linspace(0, pi, 361)
    descent    RK4 on the reduced flow, step 0.01, until |grad| < 1e-6 or 20000 steps
    endpoint   clustered at 0.02 rad; a local minimum is a cluster
    volume     each grid cell weighted by sin^6(theta) sin^5(psi) (S^7 measure
               for the two free angles); f(w) = weighted share of cells whose
               descent ends at the global minimum
    8D check   energy_and_grad(x) equals E(theta(x), psi(x)) to 1e-9 on 200
               random sphere points; analytic gradient vs central FD rel err
               < 1e-5 (as C04)

    PREMISE HOLDS at width w iff
      P1  exactly one global minimum, at angle < 0.05 from t_m, lower than
          every other local minimum by >= 0.2
      P2  f(w) <= 0.02
      P3  every grid cell whose descent ends at the global minimum has
          psi < 3 w  (the basin lies inside the declared channel sector)
      P4  the cap outside the sector is not a basin of the global minimum:
          0 cells with theta < theta_w and psi > 3 w descend to it (implied
          by P3, reported separately because it is the C20 failure)

A width failing the premise is EXCLUDED from C22's cells and listed; if
every width fails, the landscape is NOT_IN_CLASS and C22 is not run.
Nothing is retuned after the check without a new pre-registration.

## Claim C22 (registered with `scientific_method add`, UNTESTED, after the
## premise check passes and before any run)

    "On a landscape whose lowest basin is reachable only through a thin
     channel (verified exhaustively), GAS finds it at a rate different from
     uniform random at equal landscape-evaluation budget."
    experiment id: verified_channel_vs_random

## Budget: matched on landscape-function CALLS this time (budget audit)

    GAS    GASParams(max_iters = M), M in {100, 400}, default otherwise;
           rng default_rng(seed), x_init uniform from that rng; every call
           to the landscape function is counted (both term copies)
    random n_calls(GAS run, same seed) uniform sphere points, seed 1000 + seed,
           one call each — the random arm receives exactly what GAS spent
    seeds  32 per cell; cells = admitted widths x 2 budgets
    hit    angle(best_x, t_m) < r_h, r_h = 0.15 (GAS: state.best_x;
           random: argmin E over its points)

## Decision rule

Per cell two-sided Fisher exact on hits out of 32; C22 SUPPORTED iff some
cell has p < 0.05 / n_cells; direction reported per cell.

## Predictions (fixed now)

    Q1  random hit rate per cell within 2 binomial sd of 1 - (1 - f(w))^n_calls
        (the volume rate, computed from the premise check before the run)
    Q2  at the admitted thin widths (w <= 0.4) GAS < random or both zero;
        GAS does not exceed random at any admitted width at p < alpha
        (the C20/C21 direction: relaxation and thin structures do not meet)
    Q3  at w = 0.1 both arms are 0/32 at both budgets (f too small for either)

Scope: stochastic (fixed seeds), one landscape class with constants
declared above, continuous S^7; C11, C20, C21 untouched.
