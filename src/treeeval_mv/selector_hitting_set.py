"""Finite hitting-set analysis for HPR's four-term selector polynomials."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations


def selector_support(prime: int, g1: int, g2: int) -> frozenset[int]:
    """Return the nonzero coefficient positions of HPR Lemma 3.7's polynomial."""

    coefficients: dict[int, int] = {}
    terms = (
        (1, g1 * g2),
        (-1, (g1 + 1) * g2),
        (-1, g1 * (g2 + 1)),
        (1, (g1 + 1) * (g2 + 1)),
    )
    for sign, exponent in terms:
        reduced = exponent % prime
        coefficients[reduced] = coefficients.get(reduced, 0) + sign
    return frozenset(
        exponent for exponent, coefficient in coefficients.items() if coefficient
    )


@dataclass(frozen=True)
class SelectorHittingSetReport:
    prime: int
    minimum_size: int
    witness: tuple[int, ...]
    base_pairs: int
    maximum_single_exponent_coverage: int
    counting_lower_bound: int


def minimum_selector_hitting_set(prime: int) -> SelectorHittingSetReport:
    """Exhaust the minimum universal menu of coefficient positions.

    This is intended only for small-prime calibration.  The general counting
    argument needs no exhaustive search: one exponent can occur only when it
    equals one of four bilinear expressions, covering ``O(p)`` base pairs, so
    every menu hitting all ``p^2`` pairs has size ``Omega(p)``.
    """

    if prime < 3:
        raise ValueError("the HPR selector requires an odd prime")
    supports = tuple(
        selector_support(prime, g1, g2)
        for g1 in range(prime)
        for g2 in range(prime)
    )
    if any(not support for support in supports):
        raise AssertionError("HPR Lemma 3.7 predicts every support is nonempty")
    coverage = tuple(
        sum(exponent in support for support in supports) for exponent in range(prime)
    )
    maximum_coverage = max(coverage)
    counting_lower_bound = (prime * prime + maximum_coverage - 1) // maximum_coverage
    for size in range(counting_lower_bound, prime + 1):
        for candidate in combinations(range(prime), size):
            menu = frozenset(candidate)
            if all(menu & support for support in supports):
                return SelectorHittingSetReport(
                    prime=prime,
                    minimum_size=size,
                    witness=candidate,
                    base_pairs=prime * prime,
                    maximum_single_exponent_coverage=maximum_coverage,
                    counting_lower_bound=counting_lower_bound,
                )
    raise AssertionError("the full exponent set must be a hitting set")
