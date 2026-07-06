"""Partition functions of the four gases and the paper's critical constants.

Implements Theorem 4.1 (functional equations), Theorem 4.3 (the Hagedorn
ladder), Remark 4.5 (the classical Hagedorn-Frautschi-Nahm bootstrap on the
primon spectrum), and the analytic checks of Appendix A.2.

Requires ``mpmath``.
"""

from __future__ import annotations

from typing import Dict

import mpmath as mp

#: Reference values printed in the paper (Appendix A.2), used by ``verify``.
REFERENCE: Dict[str, str] = {
    "sigma_tree": "2.59738512713467162",
    "sigma_bush": "1.98847685180090810",
    "sigma_word": "1.39943332872633032",
    "kalmar_rho": "1.72864723899818362",
    "beta_boltzmann": "2.14875116590273",
    "beta_bose": "2.3072314551013476",
    "prime_zeta_2": "0.45224742004106550",
}


def prime_zeta(s):
    """The single-primon partition function P(s) = sum_p p^(-s)."""
    return mp.primezeta(s)


def z_tree(s):
    """Planar fireball gas: (1 - sqrt(1 - 4 P(s))) / 2, eq. (4.1)."""
    return (1 - mp.sqrt(1 - 4 * mp.primezeta(s))) / 2


def z_word(s):
    """Distinguishable-particle gas: P / (1 - P), eq. (4.2)."""
    p = mp.primezeta(s)
    return p / (1 - p)


def z_bush(s, ceiling: int = 60):
    """Non-planar fireball gas by the Otter cascade, eq. (4.3).

    Descends Z(s) = 1 - sqrt(1 - 2 P(s) - Z(2s)) from deep in the convergent
    region, where Z is well approximated by P.
    """
    s = mp.mpf(s)
    ladder = [s]
    while ladder[-1] < ceiling:
        ladder.append(2 * ladder[-1])
    z = mp.primezeta(ladder[-1])
    for rung in reversed(ladder[:-1]):
        z = 1 - mp.sqrt(1 - 2 * mp.primezeta(rung) - z)
    return z


def critical_exponents(dps: int = 25) -> Dict[str, object]:
    """The Hagedorn ladder of Theorem 4.3, with derivative data (Remark 4.4)."""
    with mp.workdps(dps):
        P = mp.primezeta
        sigma_tree = mp.findroot(lambda s: P(s) - mp.mpf(1) / 4, mp.mpf("2.6"))
        sigma_word = mp.findroot(lambda s: P(s) - 1, mp.mpf("1.4"))
        sigma_bush = mp.findroot(
            lambda s: 1 - 2 * P(s) - z_bush(2 * s), mp.mpf("1.95")
        )
        kalmar_rho = mp.findroot(lambda s: mp.zeta(s) - 2, mp.mpf("1.73"))
        dP_tree = mp.diff(P, sigma_tree)
        return {
            "sigma_tree": sigma_tree,
            "sigma_bush": sigma_bush,
            "sigma_word": sigma_word,
            "kalmar_rho": kalmar_rho,
            "dP_at_sigma_tree": dP_tree,
            "cusp_amplitude": mp.sqrt(-dP_tree),
            "dP_at_sigma_word": mp.diff(P, sigma_word),
            "hagedorn_temperatures": (
                mp.mpf(1),
                1 / sigma_word,
                1 / sigma_bush,
                1 / sigma_tree,
            ),
        }


def bootstrap_criticality(dps: int = 30) -> Dict[str, object]:
    """Hagedorn points of the classical bootstrap on the primon spectrum.

    Remark 4.5: Boltzmann counting is critical where P reaches Nahm's constant
    2 log 2 - 1; Bose counting solves P(s) + 2 sum_{k>=2} G(ks)/k = 2 log 2 - 1
    with G the small branch of G = P + (exp(sum_k G(ks)/k) - 1 - G).
    """
    with mp.workdps(dps):
        P = mp.primezeta
        target = 2 * mp.log(2) - 1
        beta_boltzmann = mp.findroot(lambda s: P(s) - target, mp.mpf("2.15"))

        memo: Dict[str, object] = {}

        def G(s):
            s = mp.mpf(s)
            key = mp.nstr(s, 30)
            if key in memo:
                return memo[key]
            if s > 120:
                memo[key] = P(s)
                return memo[key]
            a = tail(s)
            value_p = P(s)
            memo[key] = mp.findroot(
                lambda g: value_p + mp.e ** (g + a) - 1 - 2 * g, value_p
            )
            return memo[key]

        def tail(s):
            total, k = mp.mpf(0), 2
            while k * s <= 400:
                term = G(k * s) / k
                total += term
                if term < mp.mpf("1e-32"):
                    break
                k += 1
            return total

        beta_bose = mp.findroot(
            lambda s: P(s) + 2 * tail(s) - target, mp.mpf("2.31")
        )
        return {
            "phi_critical": target,
            "beta_boltzmann": beta_boltzmann,
            "beta_bose": beta_bose,
            "G_at_bose_criticality": (1 + P(beta_bose)) / 2,
        }


def verify(limit: int = 20000, dps: int = 25) -> Dict[str, object]:
    """Reproduce the analytic checks of Appendix A.2.

    Compares partial sums of the enumerated degeneracies against the closed and
    functional forms of Theorem 4.1 (tolerances set by the expected truncation
    tails), and every constant against the digits printed in the paper.
    """
    from primeleaf.degeneracy import brute_force_tables

    gt, gb, _ = brute_force_tables(limit)
    checks: Dict[str, bool] = {}
    with mp.workdps(dps):
        for s in (3, 4, 6):
            partial = mp.fsum(
                gt[n] * mp.mpf(n) ** (-s) for n in range(2, limit + 1)
            )
            tolerance = 10 * mp.mpf(limit) ** (mp.mpf("2.598") - s)
            checks[f"Z_tree partial sum vs closed form, s={s}"] = (
                abs(partial - z_tree(s)) < tolerance
            )
        for s in (3, 4):
            partial = mp.fsum(
                gb[n] * mp.mpf(n) ** (-s) for n in range(2, limit + 1)
            )
            tolerance = 10 * mp.mpf(limit) ** (mp.mpf("1.989") - s)
            checks[f"Z_bush partial sum vs Otter cascade, s={s}"] = (
                abs(partial - z_bush(s)) < tolerance
            )
        exponents = critical_exponents(dps)
        for name in ("sigma_tree", "sigma_bush", "sigma_word", "kalmar_rho"):
            checks[f"{name} matches the paper"] = (
                abs(exponents[name] - mp.mpf(REFERENCE[name])) < mp.mpf("1e-12")
            )
        checks["Z_bush -> 1 at its critical point"] = (
            abs(z_bush(exponents["sigma_bush"] + mp.mpf("1e-12")) - 1)
            < mp.mpf("1e-4")
        )
        checks["P(2) matches the paper"] = (
            abs(mp.primezeta(2) - mp.mpf(REFERENCE["prime_zeta_2"]))
            < mp.mpf("1e-15")
        )
        bootstrap = bootstrap_criticality(max(dps, 30))
        checks["beta_boltzmann matches the paper"] = (
            abs(bootstrap["beta_boltzmann"] - mp.mpf(REFERENCE["beta_boltzmann"]))
            < mp.mpf("1e-10")
        )
        checks["beta_bose matches the paper"] = (
            abs(bootstrap["beta_bose"] - mp.mpf(REFERENCE["beta_bose"]))
            < mp.mpf("1e-10")
        )
    return {"limit": limit, "dps": dps, "ok": all(checks.values()), "checks": checks}
