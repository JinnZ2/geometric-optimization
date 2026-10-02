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
