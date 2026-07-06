from primeleaf.primes import (
    divisors,
    factorize,
    is_prime,
    primes_up_to,
    smallest_prime_factor_sieve,
)


def test_factorize():
    assert factorize(1) == {}
    assert factorize(2) == {2: 1}
    assert factorize(360) == {2: 3, 3: 2, 5: 1}
    assert factorize(97) == {97: 1}


def test_factorize_with_sieve_agrees():
    spf = smallest_prime_factor_sieve(1000)
    for n in range(1, 1001):
        assert factorize(n, spf) == factorize(n)


def test_divisors():
    assert divisors(1) == [1]
    assert divisors(12) == [1, 2, 3, 4, 6, 12]
    assert divisors(97) == [1, 97]


def test_primes_up_to():
    assert primes_up_to(30) == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    assert primes_up_to(1) == []


def test_is_prime():
    known = set(primes_up_to(200))
    for n in range(200 + 1):
        assert is_prime(n) == (n in known)
