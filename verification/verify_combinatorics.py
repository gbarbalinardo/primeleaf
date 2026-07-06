"""Brute-force verification of the combinatorial claims in the paper (Appendix A.1).

Checks, for all n up to LIMIT:
  1. g_T(n) = Catalan(Omega(n)-1) * Omega(n)! / prod(alpha_p!)   (planar trees, OEIS A144757)
  2. g_W(n) = Omega(n)! / prod(alpha_p!)                          (words, OEIS A008480)
  3. g_B(n) unordered recursion; on n = 2^k it must reproduce the
     Wedderburn-Etherington numbers (OEIS A001190)
  4. divergence of g_B from OEIS A375120 (different diagonal convention) at n = 144

Standard library only:  python3 verify_combinatorics.py
Writes g_values.txt (n, g_T, g_B per line) for use by analytic.py.
"""
import math, sys, pathlib
from functools import lru_cache

LIMIT = 20000

# smallest-prime-factor sieve
spf = list(range(LIMIT + 1))
for i in range(2, int(LIMIT ** 0.5) + 1):
    if spf[i] == i:
        for j in range(i * i, LIMIT + 1, i):
            if spf[j] == j:
                spf[j] = i

def factor(n):
    f = {}
    while n > 1:
        p = spf[n]
        f[p] = f.get(p, 0) + 1
        n //= p
    return f

def is_prime(n):
    return n >= 2 and spf[n] == n

def divisors(n):
    ds = [1]
    for p, a in factor(n).items():
        ds = [d * p ** e for d in ds for e in range(a + 1)]
    return sorted(ds)

sys.setrecursionlimit(100000)

@lru_cache(maxsize=None)
def g_tree(n):
    """Planar binary trees (ordered + parenthesized) with shadow n."""
    t = 1 if is_prime(n) else 0
    for d in divisors(n):
        if 1 < d < n:
            t += g_tree(d) * g_tree(n // d)
    return t

@lru_cache(maxsize=None)
def g_bush(n):
    """Unordered binary trees (commutative, non-associative) with shadow n."""
    t = 1 if is_prime(n) else 0
    for d in divisors(n):
        if 1 < d < n and d * d <= n:
            if d * d == n:
                ud = g_bush(d)
                t += ud * (ud + 1) // 2      # unordered pair with repetition
            else:
                t += g_bush(d) * g_bush(n // d)
    return t

@lru_cache(maxsize=None)
def g_word(n):
    """Words in primes with product n (ordered prime factorizations)."""
    if n == 1:
        return 1  # empty word, used only inside the recursion
    return sum(g_word(n // p) for p in factor(n))

@lru_cache(maxsize=None)
def a375120(n):
    """OEIS A375120's convention: ordered pair on the diagonal split (d, d)."""
    if is_prime(n):
        return 1
    return sum(a375120(d) * a375120(n // d)
               for d in divisors(n) if 2 <= d <= n // d)

def catalan(k):
    return math.comb(2 * k, k) // (k + 1)

def multinomial(n):
    f = factor(n)
    om = sum(f.values())
    m = math.factorial(om)
    for a in f.values():
        m //= math.factorial(a)
    return m

def g_tree_formula(n):
    f = factor(n)
    return catalan(sum(f.values()) - 1) * multinomial(n)

bad_g = [n for n in range(2, LIMIT + 1) if g_tree(n) != g_tree_formula(n)]
bad_w = [n for n in range(2, LIMIT + 1) if g_word(n) != multinomial(n)]
print(f"g_T formula vs brute force, n=2..{LIMIT}: "
      f"{'ALL MATCH' if not bad_g else f'MISMATCH at {bad_g[:10]}'}")
print(f"g_W formula vs brute force, n=2..{LIMIT}: "
      f"{'ALL MATCH' if not bad_w else f'MISMATCH at {bad_w[:10]}'}")

WE = [1, 1, 1, 2, 3, 6, 11, 23, 46, 98, 207, 451, 983, 2179]  # A001190, k = 1..14
ok = all(g_bush(2 ** k) == WE[k - 1] for k in range(1, 15) if 2 ** k <= LIMIT)
print(f"g_B(2^k) vs Wedderburn-Etherington, k=1..14: {'ALL MATCH' if ok else 'MISMATCH'}")
print(f"hand checks: g_B(12)={g_bush(12)} (2), g_B(24)={g_bush(24)} (4), g_B(30)={g_bush(30)} (3)")
print(f"A375120 divergence: n=144 -> g_B={g_bush(144)} vs A375120={a375120(144)}  "
      f"(expected 41 vs 42)")

out = pathlib.Path(__file__).parent / "g_values.txt"
with open(out, "w") as fh:
    for n in range(2, LIMIT + 1):
        fh.write(f"{n} {g_tree(n)} {g_bush(n)}\n")
print(f"wrote {out}")
