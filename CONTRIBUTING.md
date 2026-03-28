# Contributing to Geometric Optimization

Thank you for your interest in contributing to this research framework.

## How to Contribute

### Bug Reports and Fixes

Open an issue on GitHub describing:
- What you expected to happen
- What actually happened
- Steps to reproduce

For fixes, submit a pull request referencing the issue.

### New Energy Terms

The framework is designed for extensibility. To add a new energy term:

1. Subclass `EnergyTerm` in `gas/energy_terms.py`
2. Implement the `compute(self, x, neighbors) -> float` method
3. Add tests in `tests/test_energy_terms.py`
4. Add a fieldlink entry in `bridges/rosetta-fieldlink.json` if the term
   corresponds to a Rosetta shape

```python
class MyEnergy(EnergyTerm):
    """Description of what geometric property this measures."""

    def compute(self, x, neighbors):
        # Your computation here
        return energy_value
```

### Performance Improvements

NumPy/SciPy optimizations, JAX ports, or algorithmic improvements are
welcome. Please include benchmark comparisons.

### Validation Studies

If you apply G-Opt to a real problem domain:
- Document the mapping from your problem to E8 space
- Report convergence behavior and solution quality
- Compare with conventional optimization where applicable

### Theoretical Extensions

Extensions to the mathematical framework (new convergence proofs,
alternative lattice structures, etc.) are welcome. Please include
derivations in a document alongside any code changes.

## Development Setup

```bash
git clone https://github.com/JinnZ2/geometric-optimization.git
cd geometric-optimization
pip install -e ".[dev]"
```

## Code Standards

- Format with `black` (88 char line length)
- Lint with `flake8`
- Type check with `mypy`
- Test with `pytest`

```bash
black .
flake8 gas/ meta_layer/ tests/
mypy gas/ meta_layer/
pytest --cov=gas --cov=meta_layer
```

All checks must pass before merging.

## Attribution

See [Contributors.md](Contributors.md) for the full list of contributors
and the collaborative research philosophy behind this project.

We contribute anonymously because the ideas matter more than the authors.
