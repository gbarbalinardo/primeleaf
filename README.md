# Prime-Leaf Tree Theory

[![CI](https://github.com/gbarbalinardo/primeleaf/actions/workflows/ci.yml/badge.svg)](https://github.com/gbarbalinardo/primeleaf/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21252898.svg)](https://doi.org/10.5281/zenodo.21252898)

Arithmetic in which multiplication does not return a number. The product of two
objects is a multiplication tree over prime leaves, and the ordinary integer is the
shadow of that tree. Read as statistical mechanics: the primes are the modes of a gas
(mode `p` has energy `log p`), a tree is one microstate (a recorded combination history),
the Riemann zeta function is the free-boson partition function, `Z = P + Z^2` is a
statistical bootstrap equation, critical exponents are inverse Hagedorn temperatures,
and ordinary arithmetic is what survives total coarse-graining. The Fundamental
Theorem of Arithmetic is exactly the statement that the coarse-grained description is
faithful.

**Paper:** [`paper/prime-leaf-tree-theory.pdf`](paper/prime-leaf-tree-theory.pdf)

## Layout

```
paper/            the paper (LaTeX + PDF), its figures, and notes
  figures/        generated figures (make with scripts/make_figures.py)
  notes/          the original informal recap and the earlier math-first draft
primeleaf/        Python package: every computation in the paper's Appendix A
scripts/          make_figures.py (regenerates the paper's figures)
tests/            fast pytest suite (the full census is marked slow)
data/             computed artifacts (tree zeros to t = 5000, the density curve)
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
| `python -m primeleaf census --t-max 5000 --sigma-min 1.25 --prime-limit 200000 --threshold 0.09 --out data/zeros_t5000.txt` | the 595 tree zeros to height 5000 | ~1 h |
| `python -m primeleaf density --sigma 1.30 1.10 1.05 1.02 1.01 1.00 --out data/density_curve.txt` | the density curve D(sigma) and total density D(1) | ~5 min |
| `python -m primeleaf dominance` | Monte Carlo q_p and the predicted zero density | ~1 min |
| `python scripts/make_figures.py` | the paper's four figures into paper/figures/ | ~1 min |
| `pytest -m slow` | full census as a test | ~10 min |

## Key results

| gas | states | inverse limiting temperature |
|---|---|---|
| multisets | classical zeta | 1 (pole) |
| words | ordered prime factorizations | 1.3994333287263303 (P = 1) |
| bushes | unordered binary trees | 1.9884768518009081 (Otter condition) |
| trees | planar binary trees | 2.5973851271346716 (P = 1/4), Z = 1/2 there |

Tree zeros (solutions of P(s) = 1/4, the complex-temperature zeros of the tree gas):
a census finds 595 of them with 0 < t < 5000, all left of sigma* = 2.597, real parts
wandering (not on a vertical line). Their density in each substrip of Re s > 1 is
given exactly by the mode-dominance formula (Theorem 5.3, proved in full), with total
density D(1) = 0.1436 (4-seed Monte Carlo cross-check). The spacing statistics show a
spectrum far more rigid than either Poisson or GUE: number variance saturates at
Sigma^2(10) = 0.67 vs 10 for Poisson, and spacings repel (min normalized gap 0.79).
This is the fingerprint of the zeros of an almost periodic function.

The classical Hagedorn-Frautschi-Nahm bootstrap transplanted onto the prime
spectrum has beta_H = 2.14875116590273 (Boltzmann counting) and
beta_H = 2.30723145510135 (Bose counting).

The unordered-tree degeneracy sequence (1, 1, 1, ..., 2 at n = 12, ..., 41 at
n = 144) is OEIS entry [A375120](https://oeis.org/A375120). Its original formula and
program counted the equal-factor split as an ordered pair (42 at n = 144 instead of
41); the verification here exposed the error, and the entry was corrected in July
2026 with a b-file to n = 10000.

## Building the paper

```sh
make paper        # or: cd paper && latexmk -pdf prime-leaf-tree-theory.tex
```

## Citing

See [`CITATION.cff`](CITATION.cff) (GitHub's "Cite this repository" button uses it).
The repository is archived on Zenodo: DOI
[10.5281/zenodo.21252898](https://doi.org/10.5281/zenodo.21252898) (all versions).
Until the paper has a journal or arXiv identifier, cite the manuscript PDF via that
DOI.
