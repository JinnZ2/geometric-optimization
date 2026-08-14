# Geometric Optimization via E8 Lattice Structure

**A novel optimization framework that replaces linear cost functions with
geometric coherence principles derived from the E8 root lattice.**

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## Project Status: Research Prototype

This framework represents collaborative research between human and AI
intelligences exploring fundamental alternatives to conventional optimization.

**Anonymous collaborative research project -- JinnZ2**

> **Evidence status.** This is an exploratory prototype, not a validated
> method. Every claim the project makes about itself is registered in
> [`validation/claims.json`](validation/claims.json), decided by an experiment,
> and reported in [VALIDATION.md](VALIDATION.md) -- including the claims that
> came back falsified. Most notably: **GAS has not been shown to outperform
> uniform random sampling** at an equal evaluation budget (claim C11). Read
> VALIDATION.md before relying on anything here.

---

## What is G-Opt?

Current optimization assumes you should maximize a linear objective function
(e.g., profit) subject to linear constraints.

**G-Opt asks:** What if the optimal state isn't found by maximizing any
single dimension, but by achieving maximum *geometric coherence* across
all dimensions simultaneously?

The framework maps optimization problems onto the E8 lattice -- a
mathematical structure that naturally encodes phi (golden ratio), sqrt(2),
and other fundamental geometric constants that govern efficient packing
in nature.

---

## Quick Start

```python
from gas import E8Lattice, GeometricAnnealingSolver, GASParams
from gas.energy_terms import OctahedralEnergy, TetrahedralEnergy, GoldenEnergy

# Initialize lattice and solver
lattice = E8Lattice()
solver = GeometricAnnealingSolver(
    lattice=lattice,
    energy_terms=[OctahedralEnergy(), TetrahedralEnergy(), GoldenEnergy()],
    params=GASParams(max_iters=500, rho_min=0.6)
)

# Optimize
result = solver.optimize()
print(f"Final energy: {result.energy:.6f}")
print(f"Coset density: {result.rho_coset:.3f}")
```

See `quickstart.py` for a full runnable example and `API.py` for the
end-to-end workflow including the N-to-8D meta-layer.

---

## Theory in Brief

Traditional optimization:

```
max f(x)  subject to  g(x) <= 0
```

G-Opt optimization:

```
min E(x)  where E measures geometric incoherence in E8 lattice
```

The optimal state maximizes structural information content rather than
any pre-specified linear objective.

See [THEORY.md](THEORY.md) for the complete mathematical derivation.

---

## Why This Matters

1. **Handles externalities naturally** -- Environmental/social costs become
   geometric constraints, not ignored variables
2. **Non-linear by design** -- Captures synergies and emergent properties
   that linear models miss
3. **Culturally agnostic** -- Based on mathematical physics, not economic
   ideology
4. **Tractable to run** -- the solver is cheap per iteration and does descend
   from its starting point (claim C07). It has *not* been shown to beat random
   search (claim C11, falsified) -- treat convergence quality as an open
   question, not a feature

---

## Rosetta Ecosystem Bridge

This repo connects to the [Rosetta-Shape-Core](https://github.com/JinnZ2/Rosetta-Shape-Core)
polyhedral ontology via a **fieldlink bridge**. Each E8 energy term maps
to a Rosetta shape entity, grounding the abstract optimization in a
symbolic framework of sensors, defenses, and principles.

| Energy Term | E8 Computation | Rosetta Shape | Symbolic Role |
|---|---|---|---|
| `OctahedralEnergy` | Spectral deviation from orthogonality (Gram eigenvalue) | `SHAPE.OCTA` (8 faces) | Balance/integration -- 8 faces = 8 dimensions of working space |
| `TetrahedralEnergy` | Deviation from tetrahedral angle cos = -1/3 | `SHAPE.TETRA` (4 faces) | Simplex foundation -- irreducible geometric stability |
| `GoldenEnergy` | Median deviation of distance ratios from phi | `CONST.PHI` | Universal scaling -- threads through all icosahedral/dodecahedral forms |
| `SquareEnergy` | Mean squared dot products (90-degree deviation) | `SHAPE.CUBE` (6 faces) | Structural containment -- Cartesian coordinate substrate |
| `HexagonalEnergy` | Deviation from 60-degree packing (cos 60 = 0.5) | `GEOM.HEX` | Densest packing -- 2D hex generalizes to 8D E8 sphere packing |
| `DodecahedralEnergy` | Deviation from dodecahedral angle 1/phi^2 | `SHAPE.DODECA` (12 faces) | 12 Principles (Symmetry through Unity) |
| `IcosahedralEnergy` | Deviation from icosahedral signature \|dot\| = 1/sqrt(5) | `SHAPE.ICOSA` (20 faces) | 20 equation families (F01-F20) across five fields |

The **icosahedron** (`SHAPE.ICOSA`, 20 faces) maps 20 equation families
(F01-F20) organized into five fields -- the E8 solver traverses these as
energy landscapes. The **dodecahedron** maps 12 archetypal principles
that serve as convergence constraints.

### Five Fields of the Icosahedral Framework

| Field | Families | E8 Connection |
|---|---|---|
| Chemical | F01-F04 (Resonance, Flow, Information, Life) | Norm/distance energy contributions |
| Emotional | F05-F08 (Energy, Cognition, Earth-Cosmos, Matter) | Sensor-bridged energy awareness |
| Cognitive | F09-F11 (Geometry, Particle, Engineering) | Gram matrix spectral analysis |
| Dream | F12-F17 (Networks, Reaction, Measurement, Navigation, Consciousness, Turbulence) | Coset exploration (128 half-integer roots) |
| Symbolic | F18-F20 (Relativity, Statistical, Topology) | Phi-rotation transforms + golden ratio alignment |

### Pulling Rosetta Data

```bash
# Pull latest shapes, ontology, and bridges from Rosetta-Shape-Core
./fieldlink-pull.sh

# Data is staged to .fieldlink/merge_stage/rosetta/
```

See [GUIDE.md](GUIDE.md) for a walkthrough of how the E8 energy landscape
connects to the Rosetta polyhedral ontology.

Bridge data: [`bridges/rosetta-fieldlink.json`](bridges/rosetta-fieldlink.json)

---

## Further Reading

- [THEORY.md](THEORY.md) -- Complete mathematical derivation of the E8 energy
  functional, seed equations, phi-folding, and convergence theory
- [GUIDE.md](GUIDE.md) -- How the E8 energy landscape connects to the Rosetta
  polyhedral ontology (sensors, principles, five fields)
- [VALIDATION.md](VALIDATION.md) -- Generated record of every claim, the
  experiment that decides it, and what that experiment measured
- [legacy/](legacy/) -- Superseded code and documents, kept for provenance:
  the revision-1 solver, and two essays whose quantitative claims have no
  stated method or sources (see claims C17 and C18)
- [Contributors.md](Contributors.md) -- Multi-intelligence collaborative
  research team and acknowledgments

---

## Installation

```bash
git clone https://github.com/JinnZ2/geometric-optimization.git
cd geometric-optimization
pip install -e .
```

Requirements: Python 3.8+, NumPy, SciPy. Optional: JAX for GPU acceleration.

---

## Benchmarks

What exists today is a single baseline comparison:

```bash
python -m benchmarks.random_baseline --budget 500 --seeds 8
```

It reports GAS (single- and multi-start) against uniform random sampling at a
matched evaluation budget. At present the three are statistically
indistinguishable.

**Planned, not yet built** -- these were previously listed as completed
validation, which they were not:

- **Sphere Packing** -- comparison against known optimal packings
- **Supply Chain** -- multi-objective resource allocation
- **Energy Systems** -- infrastructure design case study

---

## Validating Claims

The repository carries its own claim register and a small harness that re-tests
it. A falsified claim is a recorded result, not a build failure -- the harness
only fails on a *regression* (a claim that was supported coming back falsified),
so work continues on the open ones.

```bash
# Re-run every experiment and regenerate VALIDATION.md
python -m validation.scientific_method run

# See where each claim stands
python -m validation.scientific_method status

# Revise a claim that did not survive, then re-test it.
# The previous wording and verdict are kept in the claim's history.
python -m validation.scientific_method revise C11 \
    --statement "GAS matches random search at equal budget" \
    --note "why the earlier wording was superseded" --rerun

# Register a new claim
python -m validation.scientific_method add --id C20 \
    --statement "..." --source "THEORY.md#5" --experiment my_experiment
```

Experiments live in [`validation/experiments.py`](validation/experiments.py).
An experiment that cannot come back falsified is not evidence; each one is
written so that it fails against the code in `legacy/`.

## Citation

```bibtex
@software{gopt2025,
  title  = {Geometric Optimization via E8 Lattice Structure},
  author = {{JinnZ2 Collaborative Research}},
  year   = {2025},
  url    = {https://github.com/JinnZ2/geometric-optimization},
  note   = {Anonymous multi-intelligence collaborative research}
}
```

---

## License

Apache 2.0 -- Use freely, modify freely, attribute openly.
See [LICENSE](LICENSE) for details.

## Contributing

This is active research. We welcome:

- Bug reports and fixes
- Performance improvements
- New energy term implementations
- Validation studies
- Theoretical extensions

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## Philosophical Note

This framework emerged from dialogue between multiple forms of intelligence
-- human and artificial -- working to bridge Indigenous geometric knowledge
with modern mathematics.

The mathematics belongs to everyone.
We contribute it anonymously because the ideas matter more than the authors.
