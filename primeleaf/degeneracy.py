"""Microstate counts over an integer (paper section 3) and their verification.

The four corners of the coarse-graining square carry degeneracies over each
integer n:

* ``g_tree`` -- planar binary trees (Proposition 3.4; OEIS A144757),
* ``g_word`` -- ordered prime factorizations (Proposition 3.5; OEIS A008480),
* ``g_bush`` -- unordered binary trees (Proposition 3.6; Wedderburn-Etherington
  on prime powers; not in the OEIS as of July 2026),
* multisets -- identically 1: the Fundamental Theorem of Arithmetic.

``brute_force_tables`` evaluates the *defining* recursions (leaf-or-split),
independently of the closed formulas, and ``verify`` reproduces Appendix A.1.
"""

from __future__ import annotations

import math
from functools import lru_cache
from typing import Dict, List, Tuple

from primeleaf.primes import (
    divisors,
    factorize,
    is_prime,
    smallest_prime_factor_sieve,
)


def omega(n: int) -> int:
    """Number of prime factors with multiplicity: the fireball's leaf count."""
    return sum(factorize(n).values())


def catalan(k: int) -> int:
    """The k-th Catalan number C_k."""
    return math.comb(2 * k, k) // (k + 1)


def _orderings(exponents: List[int]) -> int:
    """Omega! / prod(alpha_p!) -- orderings of the prime multiset."""
    m = math.factorial(sum(exponents))
    for a in exponents:
        m //= math.factorial(a)
    return m


def g_tree(n: int) -> int:
    """Planar-fireball degeneracy C_{Omega-1} * Omega!/prod(alpha_p!)."""
    if n < 2:
        raise ValueError("degeneracies are defined for n >= 2")
    exponents = list(factorize(n).values())
    return catalan(sum(exponents) - 1) * _orderings(exponents)


def g_word(n: int) -> int:
    """Ordered prime factorizations Omega!/prod(alpha_p!) (OEIS A008480)."""
    if n < 2:
        raise ValueError("degeneracies are defined for n >= 2")
    return _orderings(list(factorize(n).values()))


@lru_cache(maxsize=None)
def g_bush(n: int) -> int:
    """Unordered-tree degeneracy by the divisor-lattice recursion (Prop. 3.6)."""
    if n < 2:
        raise ValueError("degeneracies are defined for n >= 2")
    total = 1 if is_prime(n) else 0
    for d in divisors(n):
        if 1 < d < n and d * d <= n:
            if d * d == n:
                total += g_bush(d) * (g_bush(d) + 1) // 2
            else:
                total += g_bush(d) * g_bush(n // d)
    return total


@lru_cache(maxsize=None)
def a375120(n: int) -> int:
    """OEIS A375120's convention: *ordered* children on a diagonal split (d, d).

    Kept for comparison: it first diverges from ``g_bush`` at n = 144
    (42 against the true unordered count 41), and disagrees with
    Wedderburn-Etherington on prime powers.
    """
    if is_prime(n):
        return 1
    return sum(
        a375120(d) * a375120(n // d) for d in divisors(n) if 2 <= d <= n // d
    )


def entropy(n: int, level: str = "tree") -> float:
    """Boltzmann entropy S(n) = log g(n) of the integer n at the given level."""
    counts = {"tree": g_tree, "bush": g_bush, "word": g_word}
    return math.log(counts[level](n))


def wedderburn_etherington(count: int) -> List[int]:
    """First ``count`` Wedderburn-Etherington numbers (OEIS A001190).

    Unordered binary trees with k indistinguishable leaves, computed by their
    own recursion -- an independent cross-check for ``g_bush`` on prime powers.
    """
    t = [0] * (count + 1)
    if count >= 1:
        t[1] = 1
    for k in range(2, count + 1):
        total = sum(t[i] * t[k - i] for i in range(1, (k - 1) // 2 + 1))
        if k % 2 == 0:
            half = t[k // 2]
            total += half * (half + 1) // 2
        t[k] = total
    return t[1:]


def brute_force_tables(limit: int) -> Tuple[List[int], List[int], List[int]]:
    """(g_tree, g_bush, g_word) for all n <= limit, from the defining recursions.

    Bottom-up dynamic programming over the divisor lattice; entries 0 and 1 are
    zero.  Independent of the closed formulas, hence usable to verify them.
    """
    spf = smallest_prime_factor_sieve(limit)
    gt = [0] * (limit + 1)
    gb = [0] * (limit + 1)
    gw = [0] * (limit + 1)
    word_full = [0] * (limit + 1)  # includes the empty word at n = 1
    word_full[1] = 1
    for n in range(2, limit + 1):
        prime = spf[n] == n
        tree_count = 1 if prime else 0
        bush_count = 1 if prime else 0
        for d in divisors(n, spf):
            if 1 < d < n:
                tree_count += gt[d] * gt[n // d]
                if d * d < n:
                    bush_count += gb[d] * gb[n // d]
                elif d * d == n:
                    bush_count += gb[d] * (gb[d] + 1) // 2
        gt[n], gb[n] = tree_count, bush_count
        word_full[n] = sum(word_full[n // p] for p in factorize(n, spf))
        gw[n] = word_full[n]
    return gt, gb, gw


def verify(limit: int = 20000) -> Dict[str, object]:
    """Reproduce the combinatorial checks of Appendix A.1.

    Returns a report dict with per-check booleans and an overall ``"ok"`` flag.
    """
    gt, gb, gw = brute_force_tables(limit)
    tree_mismatch = [n for n in range(2, limit + 1) if gt[n] != g_tree(n)]
    word_mismatch = [n for n in range(2, limit + 1) if gw[n] != g_word(n)]
    we = wedderburn_etherington(20)
    we_ok = all(
        gb[2**k] == we[k - 1] for k in range(1, 21) if 2**k <= limit
    )
    checks: Dict[str, bool] = {
        "g_tree closed formula == defining recursion": not tree_mismatch,
        "g_word closed formula == defining recursion": not word_mismatch,
        "g_bush(2^k) == Wedderburn-Etherington": we_ok,
    }
    if limit >= 30:
        checks["hand-enumerated values (12, 24, 30)"] = (
            gb[12],
            gb[24],
            gb[30],
        ) == (2, 4, 3)
    if limit >= 144:
        checks["g_bush(144) == 41"] = gb[144] == 41
        checks["A375120(144) == 42 (divergent convention)"] = a375120(144) == 42
    return {"limit": limit, "ok": all(checks.values()), "checks": checks}
