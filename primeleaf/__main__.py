"""Command-line entry points reproducing the paper's Appendix A.

Subcommands
-----------
verify-combinatorics   Appendix A.1: brute-force checks of the degeneracy formulas.
verify-zeta            Appendix A.2: partition-function identities and constants.
constants              Print the critical constants (Hagedorn ladder, bootstrap).
census                 Appendix A.4: tree-zero census (writes an optional zeros file).
dominance              Appendix A.4: dominance Monte Carlo and predicted density.
"""

from __future__ import annotations

import argparse
import sys
from typing import Dict, Optional, Sequence


def _print_report(report: Dict[str, object]) -> int:
    for name, ok in report["checks"].items():  # type: ignore[union-attr]
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    print("OK" if report["ok"] else "FAILED")
    return 0 if report["ok"] else 1


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="primeleaf",
        description="Reproduce the computations of the Prime-Leaf Tree Theory paper.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser(
        "verify-combinatorics",
        help="Appendix A.1: degeneracy formulas against the defining recursions",
    )
    p.add_argument("--limit", type=int, default=20000)

    p = sub.add_parser(
        "verify-zeta",
        help="Appendix A.2: partition-function identities and constants",
    )
    p.add_argument("--limit", type=int, default=20000)
    p.add_argument("--dps", type=int, default=25)

    p = sub.add_parser("constants", help="print the paper's critical constants")
    p.add_argument("--dps", type=int, default=25)

    p = sub.add_parser("census", help="Appendix A.4: tree-zero census")
    p.add_argument("--t-max", type=float, default=500.0)
    p.add_argument("--sigma-min", type=float, default=1.30)
    p.add_argument("--sigma-max", type=float, default=2.72)
    p.add_argument("--prime-limit", type=int, default=100_000)
    p.add_argument("--out", type=str, default=None, help="write zeros to this file")

    p = sub.add_parser(
        "dominance",
        help="Appendix A.4: Monte Carlo dominance probabilities and density",
    )
    p.add_argument("--sigma", type=float, default=1.30)
    p.add_argument("--samples", type=int, default=400_000)
    p.add_argument("--prime-limit", type=int, default=10_000)

    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _build_parser().parse_args(argv)

    if args.command == "verify-combinatorics":
        from primeleaf import degeneracy

        return _print_report(degeneracy.verify(args.limit))

    if args.command == "verify-zeta":
        from primeleaf import zeta

        return _print_report(zeta.verify(args.limit, args.dps))

    if args.command == "constants":
        import mpmath as mp

        from primeleaf import zeta

        exponents = zeta.critical_exponents(args.dps)
        bootstrap = zeta.bootstrap_criticality(max(args.dps, 30))
        for name in ("sigma_tree", "sigma_bush", "sigma_word", "kalmar_rho"):
            print(f"{name:16s} = {mp.nstr(exponents[name], 18)}")
        print(f"{'cusp_amplitude':16s} = {mp.nstr(exponents['cusp_amplitude'], 12)}")
        temperatures = ", ".join(
            mp.nstr(t, 10) for t in exponents["hagedorn_temperatures"]
        )
        print(f"{'hagedorn_T':16s} = {temperatures}")
        for name in ("beta_boltzmann", "beta_bose"):
            print(f"{name:16s} = {mp.nstr(bootstrap[name], 18)}")
        return 0

    if args.command == "census":
        import statistics

        import mpmath as mp

        from primeleaf import zeros as zeros_module

        result = zeros_module.census(
            args.t_max,
            (args.sigma_min, args.sigma_max),
            prime_limit=args.prime_limit,
        )
        count = len(result.zeros)
        print(f"{count} zeros with 0 < t < {args.t_max}")
        for z in result.zeros[:10]:
            print(" ", mp.nstr(z, 12))
        if count > 10:
            print(f"  ... ({count - 10} more)")
        if result.spacings:
            print(
                f"mean spacing {statistics.mean(result.spacings):.4f}"
                f"  (2*pi/log 2 = {zeros_module.LIGHTEST_MODE_PERIOD:.4f})"
            )
        if args.out:
            result.save(args.out)
            print(f"wrote {args.out}")
        return 0

    if args.command == "dominance":
        from primeleaf import zeros as zeros_module

        q, density = zeros_module.dominance_probabilities(
            args.sigma, prime_limit=args.prime_limit, samples=args.samples
        )
        first = {p: round(v, 4) for p, v in list(q.items())[:6]}
        print(f"q_p at sigma = {args.sigma}: {first}")
        print(f"predicted density: {density:.5f} zeros per unit height")
        print(f"  -> {density * 500:.1f} zeros to height 500")
        return 0

    return 2


if __name__ == "__main__":
    sys.exit(main())
