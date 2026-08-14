"""Claim registry and scientific-method harness for this repository.

See ``VALIDATION.md`` for the current record, and
``python -m validation.scientific_method --help`` for the loop that maintains it.
"""

from .experiments import Verdict, REGISTRY, experiment

__all__ = ["Verdict", "REGISTRY", "experiment"]
