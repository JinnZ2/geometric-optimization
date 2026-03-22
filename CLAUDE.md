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
├── quickstart.py           # Standalone usage example
├── API.py                  # End-to-end workflow example
├── gas/                    # Core module: Geometric Annealing Solver
│   ├── lattice.py          #   E8Lattice class (240-root system, KDTree search)
│   ├── energy_terms.py     #   Energy term ABC + concrete terms (Octahedral, Tetrahedral, Golden)
│   └── solver.py           #   GASParams, GASState, GeometricAnnealingSolver
├── meta_layer/             # N-dimensional decoder/inverse solver
│   └── decoder.py          #   ProximalGeometricDecoder (8↔N mapping via L-BFGS-B)
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
# Run tests (no test suite exists yet)
pytest

# Format code
black .

# Lint
flake8

# Type check
mypy .
```

## Architecture

### Core Flow

1. **E8Lattice** (`gas/lattice.py`) — Generates the 240-root E8 system (112 D8 roots + 128 coset vectors). Provides KDTree-based nearest-neighbor search and coset density calculation.

2. **Energy Terms** (`gas/energy_terms.py`) — Abstract `EnergyTerm` base class with implementations:
   - `OctahedralEnergy`: spectral deviation λ₁(G - I)
   - `TetrahedralEnergy`: tetrahedral angle alignment
   - `GoldenEnergy`: golden ratio (φ) alignment
   - Gradients computed via finite differences

3. **Solver** (`gas/solver.py`) — `GeometricAnnealingSolver` performs optimization with:
   - φ-rotation transforms
   - Weighted multi-term energy computation
   - Metropolis acceptance criterion
   - Adaptive annealing schedule based on coset density (ρ)
   - Convergence checking via energy tolerance

4. **Meta-layer** (`meta_layer/decoder.py`) — `ProximalGeometricDecoder` maps between N-dimensional and 8-dimensional spaces using proximal optimization (L-BFGS-B with L1 sparsity + geometric regularization).

### Key Exports (from `__init__.py`)

- `E8Lattice`, `create_test_point`
- `EnergyTerm`, `OctahedralEnergy`, `TetrahedralEnergy`, `GoldenEnergy`
- `GeometricAnnealingSolver`, `GASParams`, `GASState`

Note: Some exports declared in `__init__.py` (`SquareEnergy`, `HexagonalEnergy`, `DodecahedralEnergy`, `create_energy_suite`) are not yet implemented.

## Code Conventions

- Pure Python with NumPy/SciPy for numerics
- Dataclasses for parameter/state containers (`GASParams`, `GASState`)
- Abstract base classes for extensible components (`EnergyTerm`)
- Mathematical notation preserved in variable names (e.g., `eta_0`, `rho_0`, `tau_E`, `alpha_0`)
- 8-dimensional vectors are the native working space

## Known Gaps

- **No `__init__.py`** in `gas/` or `meta_layer/` subdirectories
- **No CLI module** — `gas/cli.py` is referenced in setup.py entry_points but does not exist
- **No test suite** — pytest is configured but no test files exist
- **No CI/CD** — no GitHub Actions or other CI configuration
- **No linter/formatter config files** — tools listed in dev deps but unconfigured
- **No pre-commit hooks**
