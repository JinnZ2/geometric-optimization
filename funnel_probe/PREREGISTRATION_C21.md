# PRE-REGISTRATION — claim C21 (geometric-optimization)

Follow-on to C20 (SUPPORTED, every significant cell GAS < random). Reading
under test (Kavik, 2026-10-02, unverified): GAS descends from where it
starts and commits to whichever basin it is in; a thin channel has almost
no volume, so descent moves away from it, while uniform sampling lands on it
at its volume rate. If GAS's avoidance carries POSITION information, sampling
where GAS does not go should find the channel better than uniform; if it is
only bad luck in descent, inverse-density sampling is at random.

Committed BEFORE C21 is added to claims.json and before any sampler code.

## C21 (to register with `validation.scientific_method add`, status UNTESTED)

    "Sampling weighted toward regions GAS avoids finds the thin channel at a
     rate different from uniform random at equal budget."
    experiment id: inverse_density_vs_random

## Build

    landscape      the C20 channel landscape, unchanged (D=1, r_b=0.5, H=2,
                   theta_w=0.8, s_w=0.1), widths w in {0.1, 0.2, 0.4, 0.8, 1.6}
    visit density  the C20 GAS arm re-run at B = 2000 with the SAME 32 seeds
                   per width, every visited x recorded (callback). Kernel
                   density on the sphere over the pooled points of that
                   width: rho(x) = sum_i exp(-|x - x_i|^2 / (2 h^2)), h = 0.3
                   (chordal distance; sphere radius sqrt 2). h is fixed now.
    sampler        rejection sampling from uniform on the sphere; a candidate
                   is accepted with probability (rho_0 / (rho(x) + rho_0)),
                   rho_0 = 0.01 * max rho over a pilot of 4096 uniform points.
                   Candidates rejected cost KDE evaluations, not energy
                   evaluations; the budget counts ENERGY evaluations only.
    budget         B = 2000 accepted samples per run, 32 seeds (rng seed
                   3000 + seed); uniform random arm = the C20 random arm
                   (seed 1000 + seed, the repo's own sphere_points)
    hits           PRIMARY  channel hit: a sample with theta in
                            [theta_w - 2 s_w, theta_w + 2 s_w] and psi < w
                            (a point in the gap of the wall)
                   SECONDARY basin hit: angle(best_x, t) < r_b (C20's criterion)
    counting       per run, 1 if any of the B samples is a hit; Fisher exact
                   two-sided per width on (inverse-density hits, uniform hits)
                   out of 32; alpha = 0.05 / 5 per criterion

## Decision rule

    C21 SUPPORTED  iff at least one width differs at p < 0.01 on the PRIMARY
                   criterion. Direction reported per width. The secondary
                   criterion is reported, not decisive.

## Cost declaration

The density is built from GAS runs costing 2000 x 32 evaluations per width.
That cost is EXCLUDED from the sampler's budget: the claim is about whether
the avoidance pattern carries position information, not about total cost.
Stated so the result is not read as an optimizer.

## Confound stated in advance

GAS never crosses the wall, so it also avoids the entire inner cap
(theta < theta_w). Inverse-density sampling therefore up-weights the basin
interior regardless of where the channel is — which is why the PRIMARY
criterion is a point in the channel itself, not the basin.

## Predictions

    R1  under the commitment reading: inverse-density > uniform on channel
        hits at w <= 0.4
    R2  under the bad-luck reading: no width differs at p < 0.01
    R3  the basin (secondary) criterion will favour inverse density at every
        width because of the stated confound; that is NOT evidence for C21
