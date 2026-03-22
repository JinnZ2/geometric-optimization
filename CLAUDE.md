# CLAUDE.md

## Project Overview

**Geometric Optimization via E8 Lattice Structure** — a Python research/scientific computing package implementing a Geometric Annealing Solver (GAS) that uses the E8 root system for optimization in 8-dimensional space, with a meta-layer for mapping between N and 8 dimensions.

- **Version:** 0.1.0 (research prototype / alpha)
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
│   ├── energy_terms.py     #   Energy term ABC + 6 concrete terms
│   ├── solver.py           #   GASParams, GASState, GeometricAnnealingSolver
│   └── cli.py              #   CLI entry point (gopt-optimize)
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
├── papers/                 # Research documents
│   └── energy.md
├── README.md
├── THEORY.md               # Full mathematical derivation (~20KB)
├── Six-Sigma.md            # Quality control analysis
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
```

## Architecture

### Core Flow

1. **E8Lattice** (`gas/lattice.py`) — Generates the 240-root E8 system (112 D8 roots + 128 coset vectors). Provides KDTree-based nearest-neighbor search and coset density calculation. `create_test_point()` generates reproducible test points on the norm-sqrt(2) sphere.

2. **Energy Terms** (`gas/energy_terms.py`) — Abstract `EnergyTerm` base class with implementations:
   - Core: `OctahedralEnergy`, `TetrahedralEnergy`, `GoldenEnergy`
   - Extended: `SquareEnergy`, `HexagonalEnergy`, `DodecahedralEnergy`
   - `create_energy_suite(include_all=False)` returns the standard or full set
   - Gradients computed via finite differences

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

49 tests across 5 test files covering lattice generation, energy terms, solver, decoder, and CLI. Run with `pytest`. Tests use `pytest.fixture` for shared lattice instances and `@pytest.mark.parametrize` for energy term variants.

## CI/CD

GitHub Actions workflow (`.github/workflows/ci.yml`) runs on push/PR to main/master:
- Matrix: Python 3.8, 3.9, 3.10, 3.11
- Steps: flake8 lint, mypy type check, pytest with coverage
