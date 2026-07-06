# Prime-Leaf Tree Theory

Arithmetic in which multiplication does not return a number. The product of two
objects is a multiplication tree over prime leaves, and the ordinary integer is the
shadow of that tree. In physics language: the primes are primons (energy log p),
trees are fireballs with recorded fusion histories, the Riemann zeta function is a
free-boson partition function, `Z = P + Z^2` is a statistical bootstrap equation,
critical exponents are inverse Hagedorn temperatures, and ordinary arithmetic is what
survives total coarse-graining. The Fundamental Theorem of Arithmetic is exactly the
statement that the coarse-grained description is faithful.

**Paper:** [`paper/prime-leaf-tree-theory.pdf`](paper/prime-leaf-tree-theory.pdf)

## Layout

```
paper/            the paper (LaTeX + PDF) and its notes
  notes/          the original informal recap and the earlier math-first draft
primeleaf/        Python package: every computation in the paper's Appendix A
tests/            fast pytest suite (the full census is marked slow)
data/             computed artifacts committed for reference (tree zeros to t = 500)
```

## Quickstart

```sh
pip install -e ".[dev]"

python -m primeleaf verify-combinatorics   # Appendix A.1  (~30 s)
python -m primeleaf verify-zeta            # Appendix A.2  (~1 min)
python -m primeleaf constants              # the critical constants
pytest                                     # fast test suite
```

## Reproducing the paper

| command | reproduces | runtime |
|---|---|---|
| `python -m primeleaf verify-combinatorics` | degeneracy formulas vs brute force, n <= 20000 | ~30 s |
| `python -m primeleaf verify-zeta` | partition-function identities, Hagedorn ladder, bootstrap constants | ~1 min |
| `python -m primeleaf census --out data/zeros_t500.txt` | the 59 tree zeros to height 500 | ~10 min |
| `python -m primeleaf dominance` | Monte Carlo q_p and the predicted zero density | ~1 min |
| `pytest -m slow` | full census as a test | ~10 min |

## Key results

| gas | states | inverse Hagedorn temperature |
|---|---|---|
| multisets | classical zeta | 1 (pole) |
| words | ordered prime factorizations | 1.3994333287263303 (P = 1) |
| bushes | unordered binary trees | 1.9884768518009081 (Otter condition) |
| trees | planar binary trees | 2.5973851271346716 (P = 1/4), Z = 1/2 there |

Tree zeros (solutions of P(s) = 1/4, the Fisher singularities of the fireball gas):
59 zeros with Re s > 1.3 and 0 < t < 500, not on a vertical line. Their density is
given exactly by the mode-dominance formula of Theorem 5.3, which predicts 58.3
zeros for that window by an independent Monte Carlo (q_2 = 0.624, q_3 = 0.148).
Eight more zeros sit below the strip, with real parts drifting down to 1.013.

The classical Hagedorn-Frautschi-Nahm bootstrap transplanted onto the primon
spectrum has beta_H = 2.14875116590273 (Boltzmann counting) and
beta_H = 2.30723145510135 (Bose counting).

The unordered-tree degeneracy sequence (1, 1, 1, ..., 2 at n = 12, ..., 41 at
n = 144) is not in the OEIS as of July 2026; A375120 is a near-miss with a different
diagonal convention. Candidate submission.

## Building the paper

```sh
make paper        # or: cd paper && latexmk -pdf prime-leaf-tree-theory.tex
```
