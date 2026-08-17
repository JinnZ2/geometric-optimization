"""Tests for the CLI module."""

import pytest

from gas.cli import main


class TestCLI:
    def test_runs_with_few_iters(self, capsys):
        main(["--max-iters", "10", "--seed", "42", "--quiet"])
        captured = capsys.readouterr()
        assert "Energy:" in captured.out

    def test_help(self):
        with pytest.raises(SystemExit) as exc:
            main(["--help"])
        assert exc.value.code == 0


def test_seed_makes_runs_reproducible(capsys):
    """--seed must reach the solver's Generator, not just the global RNG."""
    from gas.cli import main

    outputs = []
    for _ in range(2):
        main(["--max-iters", "40", "--seed", "11", "--quiet"])
        outputs.append(capsys.readouterr().out)
    assert outputs[0] == outputs[1]


def test_different_seeds_differ(capsys):
    from gas.cli import main

    main(["--max-iters", "40", "--seed", "11", "--quiet"])
    a = capsys.readouterr().out
    main(["--max-iters", "40", "--seed", "12", "--quiet"])
    assert capsys.readouterr().out != a
