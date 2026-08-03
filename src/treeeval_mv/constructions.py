"""Structured composite-modulus constructions.

The DGY generator follows Lemma 11 / HPR Appendix A.  Its low-degree weight
test is built from Lucas' theorem.  All polynomial arithmetic is exact over a
Boolean (multilinear) quotient, and the realized degree is recorded so the
degree ledger cannot silently inherit an incorrect symbolic estimate.
"""

from __future__ import annotations

from itertools import combinations, product
from math import ceil, log2, prod
from typing import Iterable, Sequence

from .model import (
    MVFamily,
    Monomial,
    Polynomial,
    PolynomialMatchingFamily,
    canonical_hpr_offdiagonal_set,
    canonical_customary_set,
    crt,
)


def coordinate_equality_hpr_family(
    *, alphabet_size: int, primes: Sequence[int], label_count: int | None = None
) -> MVFamily:
    """Return the canonical HPR family indexed by words in ``[a]^t``.

    In CRT component ``r`` the Gram entry is one exactly when the two
    labels have the same ``r``-th symbol.  Thus the diagonal CRT pattern is
    all ones, while every off-diagonal pattern contains a zero.  Each
    component equality matrix factors through ``a`` one-hot coordinates, so
    the common vector dimension is exactly ``a`` (before optional removal of
    linearly redundant coordinates in special characteristics).

    This elementary growing-prime construction is useful as a construction
    baseline: ``N=a**t`` and ``d=a``.  It is logspace uniform by direct digit
    extraction and CRT evaluation.
    """

    prime_tuple = tuple(int(prime) for prime in primes)
    if alphabet_size < 2:
        raise ValueError("alphabet_size must be at least two")
    if not prime_tuple or len(set(prime_tuple)) != len(prime_tuple):
        raise ValueError("primes must be a nonempty sequence of distinct integers")
    if any(prime < 3 for prime in prime_tuple):
        raise ValueError("HPR's selector construction requires odd prime factors")

    modulus = prod(prime_tuple)
    capacity = alphabet_size ** len(prime_tuple)
    if label_count is None:
        label_count = capacity
    if not 1 <= label_count <= capacity:
        raise ValueError("label_count must lie between one and alphabet_size**len(primes)")
    labels = tuple(product(range(alphabet_size), repeat=len(prime_tuple)))[:label_count]
    rows = []
    for label in labels:
        rows.append(
            tuple(
                crt(
                    [int(label[r] == coordinate) for r in range(len(prime_tuple))],
                    prime_tuple,
                )
                for coordinate in range(alphabet_size)
            )
        )
    ell_ceiling = ceil(log2(len(labels)))
    return MVFamily.build(
        modulus=modulus,
        u=rows,
        v=rows,
        allowed_offdiagonal=canonical_hpr_offdiagonal_set(prime_tuple),
        diagonal_value=1,
        prime_factors=prime_tuple,
        convention="hpr",
        metadata={
            "construction": "coordinate-equality CRT product",
            "alphabet_size": alphabet_size,
            "word_length": len(prime_tuple),
            "full_capacity": capacity,
            "label_bits_ceiling": ell_ceiling,
            "packed_vector_bits_ceiling": alphabet_size * ceil(log2(modulus)),
            "uniformity": "direct mixed-radix digit extraction plus CRT in log space",
        },
    )


def _clean(poly: Polynomial, modulus: int) -> Polynomial:
    return {monomial: coefficient % modulus for monomial, coefficient in poly.items() if coefficient % modulus}


def _add(left: Polynomial, right: Polynomial, modulus: int) -> Polynomial:
    result = dict(left)
    for monomial, coefficient in right.items():
        result[monomial] = (result.get(monomial, 0) + coefficient) % modulus
    return _clean(result, modulus)


def _scale(poly: Polynomial, scalar: int, modulus: int) -> Polynomial:
    return _clean({monomial: scalar * coefficient for monomial, coefficient in poly.items()}, modulus)


def _multiply(left: Polynomial, right: Polynomial, modulus: int) -> Polynomial:
    result: Polynomial = {}
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            # Variables are Boolean, so x_i^2=x_i and multiplication unions supports.
            monomial = tuple(sorted(set(left_monomial) | set(right_monomial)))
            result[monomial] = (
                result.get(monomial, 0) + left_coefficient * right_coefficient
            ) % modulus
    return _clean(result, modulus)


def _power(poly: Polynomial, exponent: int, modulus: int) -> Polynomial:
    result: Polynomial = {(): 1}
    base = poly
    power = exponent
    while power:
        if power & 1:
            result = _multiply(result, base, modulus)
        base = _multiply(base, base, modulus)
        power >>= 1
    return result


def _elementary_symmetric(variable_count: int, degree: int) -> Polynomial:
    if degree > variable_count:
        return {}
    return {tuple(indices): 1 for indices in combinations(range(variable_count), degree)}


def weight_congruence_test(
    *, variable_count: int, prime: int, exponent: int, target_weight: int
) -> Polynomial:
    """Return f with f(x)=0 iff |x|=target_weight mod prime**exponent.

    The output is valid on Boolean points and has degree at most prime**exponent-1.
    """

    if exponent < 1:
        raise ValueError("exponent must be positive")
    delta: Polynomial = {(): 1}
    for digit_index in range(exponent):
        place = prime**digit_index
        digit = _elementary_symmetric(variable_count, place)
        target_digit = (target_weight // place) % prime
        difference = _add(digit, {(): -target_digit}, prime)
        unequal = _power(difference, prime - 1, prime)
        digit_equal = _add({(): 1}, _scale(unequal, -1, prime), prime)
        delta = _multiply(delta, digit_equal, prime)
    return _add({(): 1}, _scale(delta, -1, prime), prime)


def _combine_polynomials_crt(polynomials: Sequence[Polynomial], primes: Sequence[int]) -> Polynomial:
    monomials: set[Monomial] = set()
    for polynomial in polynomials:
        monomials.update(polynomial)
    combined: Polynomial = {}
    for monomial in monomials:
        coefficient = crt([poly.get(monomial, 0) for poly in polynomials], primes)
        if coefficient:
            combined[monomial] = coefficient
    return combined


def _restrict(poly: Polynomial, retained_variables: frozenset[int]) -> Polynomial:
    return {
        monomial: coefficient
        for monomial, coefficient in poly.items()
        if set(monomial) <= retained_variables
    }


def _all_monomials(variable_count: int, maximum_degree: int) -> tuple[Monomial, ...]:
    return tuple(
        tuple(monomial)
        for degree in range(min(variable_count, maximum_degree) + 1)
        for monomial in combinations(range(variable_count), degree)
    )


def polynomial_to_mv(family: PolynomialMatchingFamily) -> MVFamily:
    u = [tuple(poly.get(monomial, 0) for monomial in family.basis) for poly in family.polynomials]
    v = [
        tuple(int(all(point[index] for index in monomial)) for monomial in family.basis)
        for point in family.points
    ]
    return MVFamily.build(
        modulus=family.modulus,
        u=u,
        v=v,
        allowed_offdiagonal=family.allowed_offdiagonal,
        diagonal_value=0,
        prime_factors=family.prime_factors,
        convention="customary",
        metadata={**family.metadata, "conversion": "DGY Lemma 36 / HPR Lemma A.4"},
    )


def dgy_grolmusz_polynomial_family(
    *, primes: Sequence[int], exponents: Sequence[int], weight: int, universe_size: int
) -> PolynomialMatchingFamily:
    """Instantiate the DGY Lemma 11 construction modeled on Grolmusz."""

    if len(primes) != len(exponents) or not primes:
        raise ValueError("primes and exponents must have equal nonzero length")
    if len(set(primes)) != len(primes):
        raise ValueError("primes must be distinct")
    if universe_size < weight:
        raise ValueError("universe_size must be at least weight")
    prime_powers = [prime**exponent for prime, exponent in zip(primes, exponents, strict=True)]
    if prod(prime_powers) <= weight:
        raise ValueError("the product of prime powers must exceed the target weight")
    component_tests = [
        weight_congruence_test(
            variable_count=universe_size,
            prime=prime,
            exponent=exponent,
            target_weight=weight,
        )
        for prime, exponent in zip(primes, exponents, strict=True)
    ]
    base = _combine_polynomials_crt(component_tests, primes)
    theoretical_degree_bound = max(prime_powers) - 1
    realized_degree = max((len(monomial) for monomial in base), default=0)
    if realized_degree > theoretical_degree_bound:
        raise AssertionError(
            f"realized degree {realized_degree} exceeds bound {theoretical_degree_bound}"
        )
    # DGY/HPR pad to all monomials through D=max p_i**e_i (one above the
    # required degree).  Preserving that padding calibrates the stated theorem.
    basis_degree = max(prime_powers)
    basis = _all_monomials(universe_size, basis_degree)
    subsets = [frozenset(values) for values in combinations(range(universe_size), weight)]
    polynomials = tuple(_restrict(base, subset) for subset in subsets)
    points = tuple(
        tuple(int(index in subset) for index in range(universe_size)) for subset in subsets
    )
    modulus = prod(primes)
    return PolynomialMatchingFamily(
        modulus=modulus,
        polynomials=polynomials,
        points=points,
        basis=basis,
        allowed_offdiagonal=canonical_customary_set(primes),
        prime_factors=tuple(primes),
        metadata={
            "construction": "DGY Lemma 11, modeled on Grolmusz; HPR Appendix A",
            "weight": weight,
            "universe_size": universe_size,
            "prime_powers": prime_powers,
            "realized_polynomial_degree": realized_degree,
            "theoretical_degree_bound": theoretical_degree_bound,
            "padded_basis_degree": basis_degree,
            "uniformity": "logspace-uniform by HPR Appendix A; this calibration materializes output",
        },
    )


def dgy_grolmusz_family(
    *, primes: Sequence[int], exponents: Sequence[int], weight: int, universe_size: int
) -> MVFamily:
    return polynomial_to_mv(
        dgy_grolmusz_polynomial_family(
            primes=primes,
            exponents=exponents,
            weight=weight,
            universe_size=universe_size,
        )
    )


def incidence_family(
    *, modulus: int, sets: Iterable[Iterable[int]], universe_size: int, allowed: Iterable[int], primes: Sequence[int]
) -> MVFamily:
    """Turn a restricted-intersection set system into an MV family."""

    normalized = [frozenset(values) for values in sets]
    rows = [tuple(int(index in values) for index in range(universe_size)) for values in normalized]
    return MVFamily.build(
        modulus=modulus,
        u=rows,
        v=rows,
        allowed_offdiagonal=allowed,
        diagonal_value=0,
        prime_factors=primes,
        convention="customary",
        metadata={"construction": "restricted-intersection incidence vectors"},
    )
