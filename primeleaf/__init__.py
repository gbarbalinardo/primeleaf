"""primeleaf: computational companion to *Prime-Leaf Tree Theory*.

Every quantitative claim in the paper — the degeneracy formulas (Propositions
3.4-3.6), the partition-function identities (Theorem 4.1), the Hagedorn ladder
(Theorem 4.3), the classical-bootstrap constants (Remark 4.5), and the tree-zero
census with its Jessen-Tornehave density check (Theorem 5.3, Appendix A) — is
reproducible from this package.

Modules
-------
``primeleaf.primes``
    Elementary prime utilities (standard library only).
``primeleaf.degeneracy``
    Microstate counts over an integer and their brute-force verification.
``primeleaf.zeta``
    Partition functions and critical constants (requires ``mpmath``).
``primeleaf.zeros``
    Tree-zero census and the dominance Monte Carlo (requires ``numpy``, ``mpmath``).

Command line
------------
``python -m primeleaf --help`` lists subcommands reproducing Appendix A.
"""

from primeleaf.degeneracy import entropy, g_bush, g_tree, g_word, omega
from primeleaf.primes import divisors, factorize, is_prime, primes_up_to

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "divisors",
    "entropy",
    "factorize",
    "g_bush",
    "g_tree",
    "g_word",
    "is_prime",
    "omega",
    "primes_up_to",
]
