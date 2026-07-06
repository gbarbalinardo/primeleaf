# Prime-Leaf Tree Theory

Arithmetic in which multiplication does not return a number: the product of two objects
is a *multiplication tree* over prime leaves, and ordinary integers are the *shadows* of
these trees. Ordinary arithmetic is recovered as the double quotient (associativity +
commutativity), and the Fundamental Theorem of Arithmetic is exactly the statement that
this collapse loses nothing numerically.

The paper is written physics-first, for readers who think in statistical mechanics:
primes are primons (energy log p), trees are fireballs with recorded fusion histories,
Z = P + Z² is a Hagedorn–Frautschi–Nahm bootstrap equation, critical exponents are
inverse Hagedorn temperatures, the complex branch points are Fisher zeros, and the
additive bridge is a process/state/expectation story. Every physics term is either an
exact structural identity or flagged as imagery (see the dictionary, Table 1).

## Files

- `prime-leaf-tree-theory.tex` / `.pdf` — the paper (working draft, physics-first).
- `notes/v1-math-first.tex` — the earlier math-first version, kept for reference.
- `notes/recap-v0.1.md` — the original informal recap the paper grew from.
- `verification/verify_combinatorics.py` — brute-force check of every counting formula
  for n ≤ 20000 (standard library only). Writes `g_values.txt`.
- `verification/analytic.py` — critical exponents, zeta-identity checks, tree-zero
  computation (`pip install mpmath numpy`; run the combinatorics script first).
- `verification/zero_census.py` — tree-zero census to height 500 plus the
  Jessen–Tornehave density check (mode-dominance Monte Carlo). Writes `zeros_t500.txt`.

## Build

```sh
latexmk -pdf prime-leaf-tree-theory.tex
```

## Key computed constants

| level | objects | critical exponent |
|---|---|---|
| multisets | classical ζ | 1 (pole) |
| words | ordered prime factorizations | 1.3994333287263303 (P = 1) |
| bushes | unordered binary trees | 1.9884768518009081 (Otter condition) |
| trees | planar binary trees | 2.5973851271346716 (P = 1/4), Z(σ*) = 1/2 |

First five "tree zeros" (solutions of P(s) = 1/4, 0 < t < 50):
1.54427 + 10.01225i, 2.14220 + 17.82247i, 2.24736 + 27.64314i,
2.15502 + 35.54071i, 2.39048 + 45.53945i — not on a vertical line.

Density theorem (§5.3, via Jessen–Tornehave): in substrips of Re s > 1 the tree zeros
have exact density (1/2π)·Σ_p log p · q_p(σ), where q_p is the probability that the
p-th mode dominates the random-phase model. Census to height 500: 59 zeros with
Re s > 1.3 vs 58.3 predicted by Monte Carlo (q_2 = 0.624, q_3 = 0.148, q_5 = 0.042);
one-per-period would give 55.2. Eight more zeros found with Re s down to 1.013 —
the real parts drift toward 1 with height (open problem 9.1).

Classical Hagedorn–Frautschi–Nahm bootstrap transplanted onto the primon spectrum
(criticality at P + corrections = 2 ln 2 − 1): β_H = 2.14875116590273 (Boltzmann
counting), β_H = 2.30723145510135 (Bose counting).

The unordered-tree degeneracy sequence (1,1,1,1,1,…,2 at n=12, …, 41 at n=144) is not
in the OEIS as of July 2026 (A375120 is a near-miss with a different diagonal
convention) — candidate submission.
