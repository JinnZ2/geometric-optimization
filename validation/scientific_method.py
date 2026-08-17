"""
A small scientific-method harness for this repository's claims.

A claim is a statement the project makes about itself. Each one is recorded in
``claims.json`` together with where it is asserted, the experiment that decides
it, and its current status. ``run`` re-executes every experiment and writes the
verdicts back; ``revise`` edits a claim that did not survive and resets it for
another attempt, keeping the prior wording in the claim's history.

A falsified claim is not a failure of the project -- it is a recorded result.
The harness only exits nonzero on a *regression*: a claim that was previously
SUPPORTED coming back FALSIFIED. Claims already known to be falsified stay
falsified without breaking the build, so work can continue on them.

Usage::

    python -m validation.scientific_method run          # re-run everything
    python -m validation.scientific_method run --claim C07
    python -m validation.scientific_method status       # print the table
    python -m validation.scientific_method revise C11 \\
        --statement "GAS matches random search at equal budget" \\
        --note "single-start GAS lost 7/8 seeds; weakened to parity"
    python -m validation.scientific_method add --id C19 --statement "..." \\
        --source "THEORY.md#5" --experiment solver_descends
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

from . import experiments

HERE = Path(__file__).parent
CLAIMS_PATH = HERE / "claims.json"
REPORT_PATH = HERE.parent / "VALIDATION.md"

SUPPORTED = "SUPPORTED"
FALSIFIED = "FALSIFIED"
UNTESTED = "UNTESTED"
UNSUPPORTED = "UNSUPPORTED"
UNFALSIFIABLE_HERE = "UNFALSIFIABLE_HERE"

BADGE = {
    SUPPORTED: "PASS",
    FALSIFIED: "FALSIFIED",
    UNTESTED: "untested",
    UNSUPPORTED: "unsupported",
    UNFALSIFIABLE_HERE: "not testable here",
}

STATUS_ORDER = [FALSIFIED, UNSUPPORTED, UNFALSIFIABLE_HERE, UNTESTED, SUPPORTED]


# --------------------------------------------------------------------------
# Registry I/O
# --------------------------------------------------------------------------
def load() -> Dict[str, Any]:
    with CLAIMS_PATH.open() as fh:
        registry: Dict[str, Any] = json.load(fh)
    return registry


def save(registry: Dict[str, Any]) -> None:
    registry["claims"].sort(key=lambda c: c["id"])
    with CLAIMS_PATH.open("w") as fh:
        json.dump(registry, fh, indent=2)
        fh.write("\n")


def find(registry: Dict[str, Any], claim_id: str) -> Dict[str, Any]:
    for c in registry["claims"]:
        if c["id"] == claim_id:
            claim: Dict[str, Any] = c
            return claim
    raise SystemExit(f"error: no claim with id {claim_id!r}")


# --------------------------------------------------------------------------
# Commands
# --------------------------------------------------------------------------
def cmd_run(args) -> int:
    registry = load()
    targets = [c for c in registry["claims"] if not args.claim or c["id"] in args.claim]
    if not targets:
        raise SystemExit("error: no matching claims")

    regressions: List[str] = []
    resolved: List[str] = []

    for claim in targets:
        fn = experiments.get(claim.get("experiment"))
        if fn is None:
            if claim.get("experiment"):
                print(
                    f"  {claim['id']}  ! experiment "
                    f"{claim['experiment']!r} is not registered"
                )
            continue

        previous = claim["status"]
        verdict = fn()
        claim["status"] = SUPPORTED if verdict.supported else FALSIFIED
        claim["evidence"] = verdict.measured
        claim["measurements"] = {k: round(v, 9) for k, v in verdict.data.items()}

        if previous == SUPPORTED and claim["status"] == FALSIFIED:
            regressions.append(claim["id"])
        if previous == FALSIFIED and claim["status"] == SUPPORTED:
            resolved.append(claim["id"])

        mark = "PASS " if verdict.supported else "FALSE"
        print(f"  {claim['id']}  [{mark}] {verdict.measured}")

    save(registry)
    write_report(registry)

    print()
    for status in STATUS_ORDER:
        n = sum(1 for c in registry["claims"] if c["status"] == status)
        if n:
            print(f"  {status:<20} {n}")

    if resolved:
        print(f"\n  resolved since last run: {', '.join(resolved)}")
    if regressions:
        print(f"\n  REGRESSION: {', '.join(regressions)} were SUPPORTED, now FALSIFIED")
        return 1
    return 0


def cmd_status(args) -> int:
    registry = load()
    for status in STATUS_ORDER:
        group = [c for c in registry["claims"] if c["status"] == status]
        if not group:
            continue
        print(f"\n{status} ({len(group)})")
        for c in group:
            print(f"  {c['id']}  {c['statement']}")
            print(f"      source: {c['source']}")
            if c.get("evidence"):
                print(f"      evidence: {c['evidence']}")
            if c.get("note"):
                print(f"      note: {c['note']}")
    return 0


def cmd_revise(args) -> int:
    """Rewrite a claim and queue it for another attempt.

    The previous wording, status and evidence are appended to the claim's
    history, so the record shows how the hypothesis changed rather than
    quietly replacing it.
    """
    registry = load()
    claim = find(registry, args.id)

    claim.setdefault("history", []).append(
        {
            "revision": claim.get("revision", 1),
            "statement": claim["statement"],
            "status": claim["status"],
            "evidence": claim.get("evidence", ""),
            "experiment": claim.get("experiment"),
            "superseded_because": args.note or "revised",
        }
    )
    claim["revision"] = claim.get("revision", 1) + 1
    if args.statement:
        claim["statement"] = args.statement
    if args.experiment:
        claim["experiment"] = args.experiment
    if args.source:
        claim["source"] = args.source
    if args.note:
        claim["note"] = args.note
    claim["status"] = UNTESTED
    claim["evidence"] = ""
    claim["measurements"] = {}

    save(registry)
    print(f"{claim['id']} -> revision {claim['revision']}, status {UNTESTED}")
    print(f"  statement: {claim['statement']}")

    if args.rerun:
        print("\nre-running:")
        args.claim = [claim["id"]]
        return cmd_run(args)
    print(
        f"\nre-run with: python -m validation.scientific_method "
        f"run --claim {args.id}"
    )
    return 0


def cmd_add(args) -> int:
    registry = load()
    if any(c["id"] == args.id for c in registry["claims"]):
        raise SystemExit(f"error: claim {args.id!r} already exists")
    registry["claims"].append(
        {
            "id": args.id,
            "statement": args.statement,
            "source": args.source,
            "kind": args.kind,
            "experiment": args.experiment,
            "status": UNTESTED,
            "evidence": "",
            "measurements": {},
            "revision": 1,
            "note": args.note or "",
            "history": [],
        }
    )
    save(registry)
    write_report(registry)
    print(f"added {args.id} (status {UNTESTED})")
    return 0


# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------
def write_report(registry: Dict[str, Any]) -> None:
    """Regenerate VALIDATION.md from the claim registry."""
    claims = sorted(
        registry["claims"],
        key=lambda c: (STATUS_ORDER.index(c["status"]), c["id"]),
    )
    lines = [
        "# Validation Record",
        "",
        "**Generated file — do not edit by hand.** Regenerate with:",
        "",
        "```bash",
        "python -m validation.scientific_method run",
        "```",
        "",
        "Every claim this project makes about itself is listed here with the",
        "experiment that decides it and what that experiment actually measured.",
        "A FALSIFIED row is a result, not a defect to hide: revise the claim with",
        "`python -m validation.scientific_method revise <ID> --statement ...` and",
        "run it again.",
        "",
        "| Status | ID | Claim | Evidence |",
        "|---|---|---|---|",
    ]
    for c in claims:
        stmt = c["statement"].replace("|", "\\|")
        ev = (c.get("evidence") or c.get("note") or "—").replace("|", "\\|")
        lines.append(f"| {BADGE[c['status']]} | {c['id']} | {stmt} | {ev} |")

    lines += ["", "## Detail", ""]
    for c in claims:
        lines.append(f"### {c['id']} — {BADGE[c['status']]}")
        lines.append("")
        lines.append(f"> {c['statement']}")
        lines.append("")
        lines.append(f"- **Asserted in:** {c['source']}")
        lines.append(f"- **Kind:** {c['kind']}")
        lines.append(f"- **Experiment:** `{c['experiment'] or 'none available'}`")
        lines.append(f"- **Revision:** {c.get('revision', 1)}")
        if c.get("evidence"):
            lines.append(f"- **Measured:** {c['evidence']}")
        if c.get("note"):
            lines.append(f"- **Note:** {c['note']}")
        if c.get("history"):
            lines.append("- **Superseded wordings:**")
            for h in c["history"]:
                lines.append(
                    f"  - r{h['revision']} ({h['status']}): "
                    f"\"{h['statement']}\" — {h['superseded_because']}"
                )
        lines.append("")

    REPORT_PATH.write_text("\n".join(lines))


# --------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="scientific_method",
        description="Run, record and revise this repository's claims.",
    )
    sub = p.add_subparsers(dest="command", required=True)

    r = sub.add_parser("run", help="run experiments and record verdicts")
    r.add_argument("--claim", action="append", help="limit to these claim ids")
    r.set_defaults(func=cmd_run)

    s = sub.add_parser("status", help="print the current claim table")
    s.set_defaults(func=cmd_status)

    v = sub.add_parser("revise", help="rewrite a claim and queue it for re-testing")
    v.add_argument("id")
    v.add_argument("--statement", help="new wording of the claim")
    v.add_argument("--experiment", help="experiment id that should decide it")
    v.add_argument("--source", help="where the claim is now asserted")
    v.add_argument("--note", help="why the previous wording was superseded")
    v.add_argument("--rerun", action="store_true", help="run it immediately")
    v.set_defaults(func=cmd_revise, claim=None)

    a = sub.add_parser("add", help="register a new claim")
    a.add_argument("--id", required=True)
    a.add_argument("--statement", required=True)
    a.add_argument("--source", required=True)
    a.add_argument("--experiment", default=None)
    a.add_argument("--kind", default="empirical")
    a.add_argument("--note", default="")
    a.set_defaults(func=cmd_add)

    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    sys.exit(main())
