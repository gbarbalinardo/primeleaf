import mpmath as mp

from primeleaf import zeta
from primeleaf.degeneracy import brute_force_tables


def test_critical_exponents_match_paper():
    exponents = zeta.critical_exponents(dps=20)
    for name in ("sigma_tree", "sigma_bush", "sigma_word", "kalmar_rho"):
        assert abs(exponents[name] - mp.mpf(zeta.REFERENCE[name])) < mp.mpf("1e-12")


def test_tree_gas_at_its_hagedorn_point():
    # Theorem 4.3: P(sigma*) = 1/4 and Z_tree(sigma*) = 1/2 exactly.
    with mp.workdps(20):
        sigma_star = mp.mpf(zeta.REFERENCE["sigma_tree"])
        assert abs(mp.primezeta(sigma_star) - mp.mpf(1) / 4) < mp.mpf("1e-15")
        assert abs(zeta.z_tree(sigma_star) - mp.mpf(1) / 2) < mp.mpf("1e-6")


def test_partial_sums_match_functional_equations():
    gt, gb, _ = brute_force_tables(2000)
    with mp.workdps(20):
        partial_tree = mp.fsum(
            gt[n] * mp.mpf(n) ** (-6) for n in range(2, 2001)
        )
        assert abs(partial_tree - zeta.z_tree(6)) < mp.mpf("1e-9")
        partial_bush = mp.fsum(
            gb[n] * mp.mpf(n) ** (-4) for n in range(2, 2001)
        )
        assert abs(partial_bush - zeta.z_bush(4)) < mp.mpf("1e-5")


def test_bush_gas_reaches_one_at_criticality():
    with mp.workdps(20):
        sigma_bush = mp.mpf(zeta.REFERENCE["sigma_bush"])
        assert abs(zeta.z_bush(sigma_bush + mp.mpf("1e-12")) - 1) < mp.mpf("1e-4")


def test_bootstrap_constants_match_paper():
    bootstrap = zeta.bootstrap_criticality(dps=25)
    assert abs(
        bootstrap["beta_boltzmann"] - mp.mpf(zeta.REFERENCE["beta_boltzmann"])
    ) < mp.mpf("1e-10")
    assert abs(
        bootstrap["beta_bose"] - mp.mpf(zeta.REFERENCE["beta_bose"])
    ) < mp.mpf("1e-9")
    # Nahm's criticality: G = (1 + P)/2 at the Bose point.
    assert abs(
        bootstrap["G_at_bose_criticality"] - mp.mpf("0.6646652572")
    ) < mp.mpf("1e-8")
