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

## C21 — inverse-density sampling (PREREGISTRATION_C21.md, commit 4316fa1; registered 8dccd85; code 7c85bf8)

Run: `python3 funnel_probe/experiment_c21.py run` from a temp directory,
which is `scientific_method run --claim C21`. Log `samples/c21_run.sample.txt`,
numbers `samples/c21_results.json`. No regression; the register reads
SUPPORTED 14, FALSIFIED 2, UNSUPPORTED 2, UNFALSIFIABLE_HERE 2, UNTESTED 1.

### What did not hold, first: all three pre-registered readings

    R1  FAILED  commitment reading (inverse-density > uniform on channel hits
                at some w <= 0.4): at w <= 0.4 both arms are at 0/32, 0/32, 0 vs 1
    R2  FAILED  bad-luck reading (no width differs): w = 0.8 differs, p = 5.9e-3
    R3  FAILED  the stated confound (basin criterion favours inverse density
                at every width): it favours UNIFORM at every width, 21/30,
                21/30, 3/30, 1/30, 1/30

So C21 is SUPPORTED by its rule -- one width differs at p < 0.01 on the
primary criterion -- and the direction is inverse-density BELOW uniform,
which is neither reading the pre-registration offered. The same shape as
C20: the rule asks "different from random" and the answer is "yes, worse".

### The table (32 seeds per arm; B = 2000 energy evaluations per run)

    w     GAS density pts   rho_0    pilot rho median   cand/accepted | channel hits inv  uni   p      | basin hits inv  uni   p
    0.1       59151          8.79        0.64              1.54        |       0/32     0/32  1.00    |     21/32    30/32  1.1e-2
    0.2       59151          8.79        0.64              1.54        |       0/32     0/32  1.00    |     21/32    30/32  1.1e-2
    0.4       62052          9.53        0.65              1.54        |       0/32     1/32  1.00    |      3/32    30/32  2.8e-12
    0.8       61290         31.5         0.44              1.28        |      22/32    31/32  5.9e-3  |      1/32    30/32  1.8e-14
    1.6       62052         55.5         0.18              1.28        |      32/32    32/32  1.00    |      1/32    30/32  1.8e-14

    w = 0.1 and 0.2 are byte-identical because GAS never reaches the channel,
    so its width never enters the energy GAS evaluates and the trajectories
    (hence the density, hence the sampler) are the same.

### What the sampler did, mechanically

The pilot median of rho (0.2-0.65) sits far below rho_0 (8.8-55), so a
uniform candidate in a region GAS never visited is accepted with
probability ~0.9: the sampler is uniform almost everywhere and carves OUT
the regions where rho exceeds rho_0 -- the GAS clusters. Two things follow.

    1. GAS clusters where its trajectories end. At w >= 0.8 GAS reaches the
       basin in 15-25 of 32 runs (C20), so the basin interior IS a GAS
       cluster and inverse density removes it: basin hits 1/32.
    2. At w <= 0.4 GAS never enters the cap (0-1 of 32), yet basin hits still
       fall 30 -> 21 and 30 -> 3. GAS stalls against the OUTER face of the
       wall, and a kernel of h = 0.3 (chordal, sphere radius 1.41) smears
       that density across a wall 0.1 wide into the cap behind it. The
       sampler cannot tell "GAS avoids the far side" from "GAS sits on the
       near side", because the kernel has no wall in it.

The channel result is the same mechanism at w = 0.8: the channel mouth is
where GAS trajectories pile up on their way in, so it is the densest GAS
region and the sampler avoids it (22/32 against uniform 31/32).

### Reading against the order's two branches

The order offered "better than random -> the avoidance carries position
information" and "at random -> bad luck in descent". The data return a
third branch: WORSE than random, because the regions GAS avoids and the
regions GAS occupies are adjacent at the scale of the kernel, and the
density records occupancy, not avoidance. GAS's avoidance of the far side
of the wall does carry position information -- the wall -- but inverting a
smoothed occupancy density does not recover it, since the information is in
the discontinuity the kernel removes. Whether a kernel narrower than the
wall (h < s_w = 0.1) would turn this around is a different sampler and a
new pre-registration; h = 0.3 was fixed before the run and is kept.

The compass line survives one more turn: inverse-density also points
south. What it is pointed at is still the wall.

Scope: stochastic (fixed seeds), one landscape class, h = 0.3, rho_0 rule
and budget as pre-registered, two implementation choices marked in the code
(pilot seed 2999; candidate batch 512, speed only); the density-building
cost (32 GAS runs x 2000 per width) excluded by declaration, so nothing
here is a statement about total cost or about C11.

## Re-scored under RULES V2 (noise band + thinness; `RULES_V2.md`, `rescore_v2.py`, nothing re-run)

Output `samples/rescore_v2.sample.txt`, numbers `samples/rescore_v2_results.json`.

    RULES V2 re-scoring -- GO C20, GAS success set per width at B=2000 (32 seeds = 32 starts, one seed each), widths [1.6, 0.8, 0.4, 0.2, 0.1]
      Rule N on GAS hits (predicted non-increasing as w falls): ['25', '15', '1', '0', '0']
         step 0: change -10  band +-7.459  ok
         step 1: change -14  band +-6.167  ok
         step 2: change -1  band +-1.918  ok
         step 3: change +0  band +-0  ok
      thinnest width with a non-empty success set: w = 0.4, share 0.031 -> T1 True ; T2 NOT_EVALUABLE (starts are random seeds, not enumerated)
      V2 OUTCOME at w = 0.4: FUNNEL_FOUND, qualifier `start-selectivity NOT_EVALUABLE`
      V2 OUTCOME at w = 0.2: NOT_FOUND_IN_RANGE (success set empty, nothing persists)
      V2 OUTCOME at w = 0.1: NOT_FOUND_IN_RANGE (success set empty, nothing persists)
      RULES_V2 prediction for GO (FUNNEL_FOUND at w = 0.4 ONLY with the qualifier; NOT_FOUND at w <= 0.2): HELD

What V2 changes here: C20's FUNNEL_FOUND narrows to w = 0.4 only (1 of 32
success, thin by T1) with T2 NOT_EVALUABLE because each start is one random
seed; at w <= 0.2 the success set is empty and nothing persists, so those
widths read NOT_FOUND_IN_RANGE. Rule N finds the width curve monotone. The
re-scoring prediction HELD. C22 (verified-thin landscape) is the
construction under which T2 could be made evaluable, by enumerating starts.

## C22 — the channel-only landscape (PREREGISTRATION_C22.md, commit 7f011fe): NOT_IN_CLASS, C22 NOT REGISTERED, GAS NOT RUN

Order: "rebuild the landscape so the basin is reachable ONLY through the
channel; verify that premise by exhaustive check BEFORE any run (new claim
C22)". Landscape `channel_landscape_v2.py` (commit 07c52f4); check
`premise_check.py`; outputs `samples/premise/`. The pre-registration's own
rule: "if every width fails, the landscape is NOT_IN_CLASS and C22 is not
run." Every width fails.

### What did not hold, first

    w     P1 unique min   P2 f(w) <= 0.02    P3 basin inside psi < 3w        P4 basin cells outside sector   PREMISE
    0.1   True (0.924)    True  (1.11e-4)    False (max psi 0.585 vs 0.300)  2767 of  7735                   FAILS
    0.2   True (0.924)    True  (3.31e-3)    False (max psi 1.169 vs 0.600)  5993                            FAILS
    0.4   True (0.924)    False (2.15e-2)    False (max psi 2.644 vs 1.200)  15191                           FAILS
    0.8   True (0.924)    False (6.54e-2)    False (max psi 3.133 vs 2.400)  7728                            FAILS
    1.6   True (0.813)    False (2.37e-1)    True  (sector is the sphere)    0                               FAILS

    ADMITTED widths: none.  landscape v2: NOT_IN_CLASS.  C22: not added to claims.json.  GAS: not run.
    361 x 361 grid, 130321 cells per width, 648 s.

The basin is NOT reachable only through the channel, and the reason is a
property of the landscape, not of the check: the bowl pushes every start
toward the antipode and the wall pushes back, and the two cancel at
theta* = 0.432 — a SHELF, a ring of equilibria just inside the wall, flat
in psi except for the channel's own exponential tail. Every start inside
the wall slides onto the shelf; the tail then pulls the shelf into the
basin from psi well beyond the sector edge (to 2x the half-width at
w = 0.1, to psi = 2.64 of pi at w = 0.4). The channel is where the shelf
DROPS into the basin, not the only way to reach it. `premise_shelf.py`
computes theta*, the psi-gradient on the shelf per width, and the drift
time from psi = pi/2 to the sector edge (3e38 time units at w = 0.1, 1.4e7
at w = 0.2, 1.4 at w = 0.4). So "reachable only through the channel" is
true at w <= 0.2 only on a timescale no optimiser runs at, and false at
w = 0.4 outright. On the GAS side this would be a random walk around the
shelf's 6-sphere until the sector cap is hit — a different experiment from
the one C22 was written for, and not one to improvise after the check.

### The check's own history (two numerical amendments, constants untouched, both before any run)

    run 1  h = 0.02 everywhere, cap 10000      78% of cells unconverged at w = 0.1, 173 "minima"
           cause: the wall's theta-curvature ~H/s_w^2 = 200 gives h*lambda = 4 > RK4's bound; cells oscillate
    run 2  h = 0.004 within 0.3 of the wall, cap 50000   killed: a 37x37 subgrid showed 73-83% unconverged at a
           CONSTANT gradient 3.00e-4 = A sin(1e-3)/2 — the theta clip at pi - 1e-3 holding cells converged
           to the antipode 1e-3 short of the pole, above TOL, forever; plus |grad| ~ 20 at w = 0.1: the
           channel's psi-curvature (D T S + H W) C / (w^2 sin^2 theta) ~ 900, h = 0.02 unstable in psi
    run 3  a cell at either clip counts as converged; h bounded by 2 / lambda_est per cell      THIS RESULT
           unconverged after 50000 steps: 497 / 1370 / 5143 / 0 / 0 cells, every one on the shelf (census)

The first reading of the run-2 stall was the shelf; the bulk was the clip.
Both are recorded in `premise_shelf.py`'s docstring in the order found.
The amendments change where the integrator steps and what counts as
settled; no landscape constant, no premise threshold and no claim text
moved. The remaining unconverged cells are shelf drifters at gradients
between TOL and 1e-3, counted beside P3 so the t -> inf limit is not read
as a pass; they are already outside the sector, so including them could
only make P3 fail harder.

### Two things recorded against this session's own work

The pre-registered Q1 compared the random arm to 1 - (1 - f(w))^n with
f(w) the basin's VOLUME fraction, while the hit criterion is a cap of
angular radius r_h = 0.15 around t_m — a different fraction, f_h =
2.44e-7 on S^7. Written into `premise_check.py` and `experiment_c22.py`
before any run (Q1 as registered beside Q1' corrected), and now moot.
And `experiment_c22.py` (commit 1f7ba1f) is code that has never executed:
it reads the admitted widths from `premise_check.json`, finds none, and is
left in place as the experiment C22 would have been.

Scope: deterministic gradient flow on the reduced (theta, psi) landscape;
the premise is a statement about gradient-flow basins, and the "only
through the channel" property GAS would experience is a Metropolis
question on the shelf that this check does not ask. Nothing here is a
GAS result; C20 and C21 stand as recorded above.
