# funnel_probe — results (geometric-optimization)

Pre-registration: commit 31d5523. C20 registered UNTESTED (commit after it,
`validation.scientific_method add`), run with
`python3 funnel_probe/experiment_c20.py run` from a temp directory, which is
`scientific_method run --claim C20` with the experiment registered at
import. Full run log: `samples/c20_run.sample.txt`; per-cell numbers:
`samples/c20_results.json`. Step 1 (noise): `NOISE_RESULTS.md`.

## What did not hold, first

Nothing pre-registered failed — all four predictions P1–P4 held — and that
is the thing to say carefully: the predictions were written expecting GAS to
lose at thin widths, and it did. What did not hold is the reading of C20
one might reach for. C20 is SUPPORTED, and the direction in every
significant cell is **GAS below random**:

    width  budget  GAS   random   Fisher p    direction
    0.1     500    0/32  11/32    3.5e-04     GAS<random
    0.1    2000    0/32  30/32    6.1e-16     GAS<random
    0.2     500    0/32  11/32    3.5e-04     GAS<random
    0.2    2000    0/32  30/32    6.1e-16     GAS<random
    0.4     500    1/32  11/32    2.7e-03     GAS<random
    0.4    2000    1/32  30/32    1.8e-14     GAS<random
    0.8     500   11/32  11/32    1.0         tie
    0.8    2000   15/32  30/32    6.8e-05     GAS<random
    1.6     500   21/32  11/32    2.4e-02     GAS>random   (not < 0.005)
    1.6    2000   25/32  30/32    1.5e-01     GAS<random   (not < 0.005)
    NO_WALL 500   32/32  11/32    6.3e-09     GAS>random   (control, not a cell)
    NO_WALL 2000  32/32  30/32    4.9e-01     GAS>random   (control)

    7 of 10 cells at p < 0.05/12 = 0.0042; all seven GAS<random.

A second thing to record against GAS rather than for it: the solver's own
convergence test (`optimize`, after step 50: energy stable AND
rho_coset > rho_min) stopped some runs at iteration 52, so the mean GAS
budget actually spent was 472–486 of 500 and 1848–1939 of 2000 while random
always spent the full budget. That is the engine's behaviour and is left in.

## Outcome: FUNNEL_FOUND (on this landscape class, for GAS's dynamics)

The set of starts from which GAS reaches the lowest basin narrows with the
channel width and does not vanish: 25/32 at w = 1.6, 15/32 at 0.8, 1/32 at
0.4, 0/32 at 0.2 and 0.1 (B = 2000). Random is flat at 30/32 because the
basin's S^7 cap fraction (~9e-4) does not depend on the wall; the
pre-registered expectation was ~0.84 at B = 2000 and the measured rate is
0.94. The no-wall control is 32/32, so GAS's gradient flow does reach the
basin when nothing is in the way: the width is what it sees.

Scope: GAS is a Metropolis annealer, so each run is a noisy trajectory with
a fixed seed; parameter set as pre-registered (D=1, r_b=0.5, H=2,
theta_w=0.8, s_w=0.1, default GASParams except max_iters); continuous
8-dimensional sphere. This is one landscape class with a hand-built wall.
It says nothing about E8 energy terms, and C11 stays exactly as recorded
(FALSIFIED, revision 1, untouched by this run).

## Claim register state after the run

C20 SUPPORTED (evidence string in `claims.json`); no regression; the
status table printed by the repo's tool: FALSIFIED 2, UNSUPPORTED 2,
UNFALSIFIABLE_HERE 2, UNTESTED 1, SUPPORTED 13.

## Budget match and the per-cell table (follow-on order, 2026-10-02)

`budget_audit.py` re-runs the C20 grid with the same seeds, counting every
call to the landscape function per arm and recording GAS trajectories.
Output: `samples/budget_audit.sample.txt`, `samples/budget_audit.json`.

### What did not hold: "equal budget" was equal in NEW POINTS, not in evaluations

    one GAS iteration  = 1 new point (the proposal)
                       + 1 re-evaluation of the current point
                       + 1 gradient at the current point
                       x 2 landscape copies (terms() supplies two)
                       = 6.0 landscape-function calls, 1 distinct new point
    one random sample  = 1 call, 1 distinct point

    B      GAS iters (mean)  GAS calls  GAS distinct pts  random calls = pts
    500        472-486        2836-2918     121-184            500
    2000      1848-1939      11093-11637   1337-1431          2000

So GAS received about 5.7x the function calls AND a gradient at every point
random never got, and still spent FEWER distinct points because the
Metropolis chain revisits and the convergence test stopped some runs early.
Under either accounting GAS is not budget-starved relative to random, so the
C20 direction (GAS < random at thin widths) is not a budget artefact. The
headline "equal budget" should be read as "equal max_iters".

### The table (Fisher two-sided; adjusted p = 12 x p, Bonferroni over 12 cells)

    width  B     GAS hits  rnd hits  Fisher p   adj p     in-channel GAS  rnd   cone fraction
    0.1    500      0/32     11/32   3.47e-04   3.47e-03        0/32    0/32     1.55e-07
    0.1    2000     0/32     30/32   6.12e-16   6.12e-15        0/32    0/32     1.55e-07
    0.2    500      0/32     11/32   3.47e-04   3.47e-03        0/32    0/32     9.75e-06
    0.2    2000     0/32     30/32   6.12e-16   6.12e-15        0/32    0/32     9.75e-06
    0.4    500      1/32     11/32   2.65e-03   2.65e-02        1/32    0/32     5.79e-04
    0.4    2000     1/32     30/32   1.79e-14   1.79e-13        1/32    1/32     5.79e-04
    0.8    500     11/32     11/32   1.00e+00   1.00e+00       12/32   21/32     2.74e-02
    0.8    2000    15/32     30/32   6.79e-05   6.79e-04       16/32   31/32     2.74e-02
    1.6    500     21/32     11/32   2.37e-02   2.37e-01       22/32   32/32     5.27e-01
    1.6    2000    25/32     30/32   1.48e-01   1.00e+00       26/32   32/32     5.27e-01

    32 seeds per arm per cell; "hits" = basin hit, angle(best_x, t) < r_b;
    "in-channel" = the run put at least one evaluated point inside the gap of
    the wall (theta within 2 s_w of theta_w, psi < w).

### Width curve against the candidate reading

The reading: random lands on the channel at its volume rate; GAS commits to
the basin it starts in and falls off faster than linear in w. The curve
does not come out that way, and the reason is informative:

    w      cone fraction   expected uniform runs touching channel (B=2000)   measured rnd   GAS
    0.1    1.55e-07        1 - (1 - 1.55e-7)^2000 = 3.1e-4  -> 0.01 of 32        0/32     0/32
    0.2    9.75e-06        0.019                            -> 0.6 of 32         0/32     0/32
    0.4    5.79e-04        0.69                             -> 22 of 32          1/32     1/32
    0.8    2.74e-02        1.00                             -> 32 of 32         31/32    16/32
    1.6    0.527           1.00                             -> 32 of 32         32/32    26/32

Random is at or below the volume rate (at w = 0.4 it is well below it,
1/32 against an expected 22 of 32; the cone fraction is the S^6 azimuthal
gap over the whole theta range, while the channel is also confined to the
2 s_w band in theta, so the true volume is smaller than the printed cone
fraction by roughly the band's share of the sphere). GAS touches the channel
at the SAME rate as random at w <= 0.4 (0/0, 0/0, 1/1) and below it at
w >= 0.8. So at thin widths neither arm sees the channel at all and the
basin-hit difference (0 vs 30 at B=2000) is produced ELSEWHERE: random
reaches the basin by sampling its ~9e-4 cap directly, from the far side of
the wall; GAS never crosses the wall and so never samples the cap. The
channel is not the route by which random wins. The commitment reading
survives in a weaker form -- GAS stays on the wall's outer side -- but
"random finds the channel at its volume rate" is not what the data show,
and "GAS falls faster than linear in w" cannot be read off a curve where
both arms are at zero below w = 0.4.

The compass line stands on the C20 table: at every thin width GAS points
away from the basin at Bonferroni. What it points at is the far side of a
wall, not the mouth of a channel.
