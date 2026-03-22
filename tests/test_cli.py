"""Tests for the CLI module."""

import pytest

from gas.cli import main


class TestCLI:
    def test_runs_with_few_iters(self, capsys):
        code = main(["--max-iters", "10", "--seed", "42", "--quiet"])
        captured = capsys.readouterr()
        assert "Energy:" in captured.out

    def test_help(self):
        with pytest.raises(SystemExit) as exc:
            main(["--help"])
        assert exc.value.code == 0
