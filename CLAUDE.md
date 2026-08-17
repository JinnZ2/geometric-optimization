# CLAUDE.md

## Project Overview

**Geometric Optimization via E8 Lattice Structure** — a Python research/scientific computing package implementing a Geometric Annealing Solver (GAS) that uses the E8 root system for optimization in 8-dimensional space, with a meta-layer for mapping between N and 8 dimensions.

- **Version:** 0.2.0 (research prototype / alpha — see VALIDATION.md)
- **License:** Apache 2.0
- **Python:** 3.8+

## Repository Structure

```
geometric-optimization/
├── __init__.py            # Package exports (top-level)
├── setup.py               # Package config (setuptools)
├── pyproject.toml          # Tool config (black, flake8, mypy, pytest)
├── quickstart.py           # Standalone usage example
├── API.py                  # End-to-end workflow example
├── gas/                    # Core module: Geometric Annealing Solver
│   ├── __init__.py         #   Subpackage exports
│   ├── lattice.py          #   E8Lattice class (240-root system, KDTree search)
│   ├── energy_terms.py     #   Energy term ABC + 7 concrete terms (revision 2)
│   ├── solver.py           #   GASParams, GASState, GeometricAnnealingSolver
│   └── cli.py              #   CLI entry point (gopt-optimize)
├── validation/             # Claim register + scientific-method harness
│   ├── claims.json         #   Every claim, its source, status and history
│   ├── experiments.py      #   Experiments that decide each claim
│   └── scientific_method.py #  CLI: run / status / revise / add
├── benchmarks/
│   └── random_baseline.py  #   GAS vs uniform random sampling
├── legacy/                 # Superseded code and documents (see legacy/README.md)
│   ├── gas/                #   Revision-1 energy_terms.py and solver.py
│   ├── Six-Sigma.md        #   Moved: claim C17, no stated method
│   └── papers/energy.md    #   Moved: claim C18, no citations
├── meta_layer/             # N-dimensional decoder/inverse solver
│   ├── __init__.py         #   Subpackage exports
│   └── decoder.py          #   ProximalGeometricDecoder (8↔N mapping via L-BFGS-B)
├── tests/                  # Test suite (pytest)
│   ├── test_lattice.py
│   ├── test_energy_terms.py
│   ├── test_solver.py
│   ├── test_decoder.py
│   └── test_cli.py
├── .github/workflows/ci.yml  # CI pipeline (lint + type check + test)
├── .pre-commit-config.yaml   # Pre-commit hooks (black, flake8, whitespace)
├── README.md
├── THEORY.md               # Full mathematical derivation (~20KB)
├── VALIDATION.md           # GENERATED - claim/evidence record, do not hand-edit
├── CITATION.cff
└── Contributors.md
```

## Build & Install

```bash
# Install in editable mode
pip install -e .

# With dev tools
pip install -e ".[dev]"

# With GPU support (JAX/CUDA)
pip install -e ".[gpu]"

# With visualization
pip install -e ".[viz]"
```

## Dependencies

| Category | Packages |
|----------|----------|
| **Core** | `numpy>=1.20.0`, `scipy>=1.7.0` |
| **Dev** | `pytest>=6.0`, `pytest-cov>=2.12`, `black>=21.0`, `flake8>=3.9`, `mypy>=0.900` |
| **GPU** | `jax[cuda]`, `jaxlib` |
| **Viz** | `matplotlib`, `seaborn`, `plotly` |

## Commands

```bash
# Run tests
pytest

# Run tests with coverage
pytest --cov=gas --cov=meta_layer --cov-report=term-missing

# Format code
black .

# Lint
flake8 gas/ meta_layer/ tests/

# Type check
mypy gas/ meta_layer/

# CLI optimizer
gopt-optimize --max-iters 500 --seed 42

# Re-decide every recorded claim and regenerate VALIDATION.md
python -m validation.scientific_method run
python -m validation.scientific_method status
python -m validation.scientific_method revise C11 --statement "..." --rerun

# Baseline benchmark
python -m benchmarks.random_baseline --budget 500 --seeds 8
```

## Scientific Integrity Workflow

This repository tracks its own claims. Before asserting anything about what the
code achieves, check `VALIDATION.md`. When changing behaviour:

1. Run `python -m validation.scientific_method run` — it fails only on a
   *regression* (a SUPPORTED claim coming back FALSIFIED).
2. If a claim's wording no longer matches what is being tested, use `revise`
   rather than editing `claims.json` by hand — it preserves the prior wording
   and the reason it was superseded.
3. A new capability claim needs a new experiment that can actually fail.
   `TestFalsification` in the test suite is written to fail against `legacy/`.

**Known open results:** GAS has not been shown to beat random search (C11).
THEORY.md 9.2 is unproved (C13). Do not describe the method as validated.

## Architecture

### Core Flow

1. **E8Lattice** (`gas/lattice.py`) — Generates the 240-root E8 system (112 D8 roots + 128 coset vectors). Provides KDTree-based nearest-neighbor search and coset density calculation. `create_test_point()` generates reproducible test points on the norm-sqrt(2) sphere.

2. **Energy Terms** (`gas/energy_terms.py`) — Abstract `EnergyTerm` base class with implementations:
   - Core: `OctahedralEnergy`, `TetrahedralEnergy`, `GoldenEnergy`
   - Extended: `SquareEnergy`, `HexagonalEnergy`, `DodecahedralEnergy`,
     `IcosahedralEnergy`
   - `create_energy_suite(include_all=False)` returns the standard or full set
   - Each term scores x-to-neighbour cosines against a `target_cosines` set,
     kernel-weighted by proximity; **analytic** gradients, verified against
     finite differences (claim C04)

3. **Solver** (`gas/solver.py`) — `GeometricAnnealingSolver` performs optimization with:
   - phi-rotation transforms
   - Weighted multi-term energy computation
   - Metropolis acceptance criterion
   - Adaptive annealing schedule based on coset density (rho)
   - Convergence checking via energy tolerance

4. **Meta-layer** (`meta_layer/decoder.py`) — `ProximalGeometricDecoder` maps between N-dimensional and 8-dimensional spaces using proximal optimization (L-BFGS-B with L1 sparsity + geometric regularization).

5. **CLI** (`gas/cli.py`) — Command-line interface via `gopt-optimize` entry point. Supports `--max-iters`, `--seed`, `--rho-min`, `--all-terms`, `--quiet`.

### Key Exports (from `__init__.py`)

- `E8Lattice`, `create_test_point`
- `EnergyTerm`, `OctahedralEnergy`, `TetrahedralEnergy`, `GoldenEnergy`, `SquareEnergy`, `HexagonalEnergy`, `DodecahedralEnergy`, `create_energy_suite`
- `GeometricAnnealingSolver`, `GASParams`, `GASState`

## Code Conventions

- Pure Python with NumPy/SciPy for numerics
- `from __future__ import annotations` for forward references in type hints
- Dataclasses for parameter/state containers (`GASParams`, `GASState`)
- Abstract base classes for extensible components (`EnergyTerm`)
- Mathematical notation preserved in variable names (e.g., `eta_0`, `rho_0`, `tau_E`, `alpha_0`)
- 8-dimensional vectors are the native working space
- Formatting: black (88 char line length), flake8 for linting, mypy for type checking

## Testing

105 tests across 5 test files covering lattice generation, energy terms, solver, decoder, and CLI. Run with `pytest`. Tests use `pytest.fixture` for shared lattice instances and `@pytest.mark.parametrize` for energy term variants.

`TestFalsification` classes in `test_energy_terms.py` and `test_solver.py` encode
one test per defect found in the revision-1 audit; each is written to fail
against the corresponding file in `legacy/gas/`. Do not weaken them — a
regression test that cannot fail is not evidence.

## CI/CD

GitHub Actions workflow (`.github/workflows/ci.yml`) runs on push/PR to main/master:
- Matrix: Python 3.10, 3.11, 3.12
- Steps: flake8 lint, mypy type check, pytest with coverage, claim re-test,
  and a check that `VALIDATION.md` is not stale
