"""Elementary prime-number utilities.

Standard library only, sized for the ranges used in the paper's verifications
(n up to a few times 10**5).  For bulk work over every n <= limit, build a
smallest-prime-factor sieve once and pass it to the other helpers.
"""

from __future__ import annotations

from typing import Dict, List, Optional


def smallest_prime_factor_sieve(limit: int) -> List[int]:
    """Return ``spf`` with ``spf[n]`` the smallest prime factor of ``n >= 2``."""
    if limit < 1:
        raise ValueError("limit must be >= 1")
    spf = list(range(limit + 1))
    i = 2
    while i * i <= limit:
        if spf[i] == i:
            for j in range(i * i, limit + 1, i):
                if spf[j] == j:
                    spf[j] = i
        i += 1
    return spf


def primes_up_to(limit: int) -> List[int]:
    """All primes p <= limit, ascending."""
    if limit < 2:
        return []
    spf = smallest_prime_factor_sieve(limit)
    return [n for n in range(2, limit + 1) if spf[n] == n]


def factorize(n: int, spf: Optional[List[int]] = None) -> Dict[int, int]:
    """Prime factorization ``{p: alpha_p}`` of ``n >= 1`` (empty dict for 1).

    Uses the sieve when provided and applicable, trial division otherwise.
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    factors: Dict[int, int] = {}
    if spf is not None and n < len(spf):
        while n > 1:
            p = spf[n]
            factors[p] = factors.get(p, 0) + 1
            n //= p
        return factors
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def divisors(n: int, spf: Optional[List[int]] = None) -> List[int]:
    """All positive divisors of ``n``, ascending."""
    ds = [1]
    for p, a in factorize(n, spf).items():
        ds = [d * p**e for d in ds for e in range(a + 1)]
    return sorted(ds)


def is_prime(n: int, spf: Optional[List[int]] = None) -> bool:
    """Primality by sieve lookup or trial division (fine for these ranges)."""
    if n < 2:
        return False
    if spf is not None and n < len(spf):
        return spf[n] == n
    return list(factorize(n).values()) == [1]
