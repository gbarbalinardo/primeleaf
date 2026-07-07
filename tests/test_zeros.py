import pytest

from primeleaf import zeros


def test_first_tree_zero():
    # Table 3, row 1: the census in a small box finds exactly the first zero.
    result = zeros.census(t_max=12.0, sigma_range=(1.4, 2.7), prime_limit=20_000)
    assert len(result.zeros) == 1
    z = result.zeros[0]
    assert abs(z.real - 1.54427014392) < 1e-8
    assert abs(z.imag - 10.0122461654) < 1e-8


def test_census_bookkeeping():
    result = zeros.census(t_max=30.0, sigma_range=(1.4, 2.7), prime_limit=20_000)
    # Table 3: two zeros below t = 30 (t = 10.01 and t = 17.82; third is 27.64).
    assert result.counting(25.0) == 2
    assert len(result.spacings) == len(result.zeros) - 1
    assert all(s > 0 for s in result.spacings)


def test_census_save_load_roundtrip(tmp_path):
    result = zeros.census(t_max=12.0, sigma_range=(1.4, 2.7), prime_limit=20_000)
    path = tmp_path / "zeros.txt"
    result.save(path)
    loaded = zeros.Census.load(path)
    assert len(loaded.zeros) == len(result.zeros)
    assert abs(loaded.zeros[0] - result.zeros[0]) < 1e-12


def test_dominance_probabilities_sane():
    q, density = zeros.dominance_probabilities(
        1.30, prime_limit=2000, samples=20_000, track=20, seed=7
    )
    # Theorem 5.3(iii) at sigma = 1.30: q_2 ~ 0.62, density ~ 0.117 (Appendix A.4).
    assert 0.5 < q[2] < 0.75
    assert q[2] > q[3] > q[5]
    assert 0.08 < density < 0.16


def test_density_prediction_curve_is_increasing():
    # Corollary 5.7: D(sigma) increases as sigma decreases toward 1.
    d130 = zeros.density_prediction(
        1.30, prime_limit=2000, samples=30_000, track=40, seed=3
    )
    d100 = zeros.density_prediction(
        1.00, prime_limit=2000, samples=30_000, track=40, seed=3
    )
    assert d100["density"] > d130["density"]
    # total density D(1) ~ 0.1436 (4-seed cross-check, Appendix A.4)
    assert 0.12 < d100["density"] < 0.17
    assert d100["tail"] >= 0.0


@pytest.mark.slow
def test_full_census_matches_paper():
    # Appendix A.4: 59 zeros with Re s > 1.3 and 0 < t < 500.
    result = zeros.census(t_max=500.0)
    assert len(result.zeros) == 59
