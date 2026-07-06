import math

from primeleaf.degeneracy import (
    a375120,
    brute_force_tables,
    entropy,
    g_bush,
    g_tree,
    g_word,
    omega,
    verify,
    wedderburn_etherington,
)

LIMIT = 2000


def test_closed_formulas_match_defining_recursions():
    gt, gb, gw = brute_force_tables(LIMIT)
    for n in range(2, LIMIT + 1):
        assert gt[n] == g_tree(n)
        assert gw[n] == g_word(n)
        assert gb[n] == g_bush(n)


def test_paper_examples():
    # Example 3.8: the four levels over n = 12.
    assert (g_tree(12), g_bush(12), g_word(12)) == (6, 2, 3)
    assert g_tree(24) == 20
    assert g_bush(24) == 4
    assert g_bush(30) == 3
    assert omega(12) == 3
    assert math.isclose(entropy(12), math.log(6))


def test_wedderburn_etherington_on_prime_powers():
    we = wedderburn_etherington(10)
    assert we == [1, 1, 1, 2, 3, 6, 11, 23, 46, 98]
    for k in range(1, 11):
        if 2**k <= LIMIT:
            assert g_bush(2**k) == we[k - 1]
    # The count is prime-blind (Proposition 8.1(ii)): same values on 3^k.
    for k in range(1, 7):
        assert g_bush(3**k) == we[k - 1]


def test_a375120_divergence_at_144():
    # Remark 3.7: the OEIS near-miss counts ordered children on (d, d) splits.
    assert g_bush(144) == 41
    assert a375120(144) == 42


def test_shape_blindness():
    # Proposition 8.1(ii): degeneracies see only the exponent profile.
    for a, b in [(12, 18), (12, 20), (2**3 * 3, 5**3 * 7)]:
        assert g_tree(a) == g_tree(b)
        assert g_bush(a) == g_bush(b)
        assert g_word(a) == g_word(b)


def test_verify_report():
    report = verify(500)
    assert report["ok"], report
