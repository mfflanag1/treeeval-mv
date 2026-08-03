"""Cohomological barrier for branch-aware stateless CRT selectors."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import gcd, prod
from typing import Sequence

from .model import crt


@dataclass(frozen=True)
class BranchAwareBarrierReport:
    modulus: int
    primes: tuple[int, ...]
    exact_selector_impossible: bool
    telescoping_prime: int
    nonzero_in_modulus: bool
    arbitrary_finite_shift_sets: bool


def unit_edge_potential(*, cycle_length: int, modulus: int) -> tuple[int, ...]:
    """Construct a cyclic potential whose adjacent differences are units.

    This positive calibration explains why a one-step *projective* selector
    can exist even though exact normalization cannot.  Both arguments require
    odd cycle length and odd modulus.
    """

    if cycle_length < 3 or cycle_length % 2 == 0 or modulus % 2 == 0:
        raise ValueError("cycle length and modulus must be odd, with length >= 3")
    increments = [1, 1, -2]
    for _ in range((cycle_length - 3) // 2):
        increments.extend((1, -1))
    if sum(increments) != 0 or len(increments) != cycle_length:
        raise AssertionError("invalid zero-sum unit increment sequence")
    if any(gcd(value, modulus) != 1 for value in increments):
        raise ValueError("the construction requires every increment to be a unit")
    values = [0]
    for increment in increments[:-1]:
        values.append((values[-1] + increment) % modulus)
    return tuple(values)


def verify_unit_edge_potential(values: Sequence[int], modulus: int) -> bool:
    """Check cyclic restoration and unit adjacent differences exactly."""

    if not values:
        return False
    return all(
        gcd((values[(index + 1) % len(values)] - values[index]) % modulus, modulus)
        == 1
        for index in range(len(values))
    )


def branch_aware_exact_barrier(primes: Sequence[int]) -> BranchAwareBarrierReport:
    """Return the exact torus-telescoping obstruction.

    Let a one-dimensional branch estimator use any finite set of additive
    shifts ``S`` and arbitrary shift-specific functions ``A_s`` of the current
    full CRT inner product.  Sum the target equation around a cycle in the
    first CRT component.  Reindexing separately for every shift turns the
    left side into the sum of the single-second-component-delta equation,
    hence zero.  The target right side sums to ``p_1`` in ``Z_m``.  Since
    ``m`` has another coprime factor, ``p_1`` is nonzero modulo ``m``.

    Any two-input selector restricts to this estimator after fixing the second
    input to a target edge, so the obstruction applies even when a branch
    weight sees both current inner products separately.
    """

    prime_tuple = tuple(primes)
    if len(prime_tuple) < 2 or len(set(prime_tuple)) != len(prime_tuple):
        raise ValueError("at least two distinct CRT primes are required")
    modulus = prod(prime_tuple)
    prime = prime_tuple[0]
    return BranchAwareBarrierReport(
        modulus=modulus,
        primes=prime_tuple,
        exact_selector_impossible=True,
        telescoping_prime=prime,
        nonzero_in_modulus=prime % modulus != 0,
        arbitrary_finite_shift_sets=True,
    )


def multiplicative_unit_subgroups(modulus: int) -> tuple[tuple[int, ...], ...]:
    """Enumerate every multiplicatively closed unit scale set containing one."""

    units = tuple(value for value in range(modulus) if gcd(value, modulus) == 1)
    subgroups: list[tuple[int, ...]] = []
    for size in range(1, len(units) + 1):
        for candidate in combinations(units, size):
            values = set(candidate)
            if 1 not in values:
                continue
            if not all(a * b % modulus in values for a in values for b in values):
                continue
            if not all(pow(a, -1, modulus) in values for a in values):
                continue
            subgroups.append(candidate)
    return tuple(subgroups)


def crt_binary_shifts(primes: Sequence[int]) -> tuple[int, ...]:
    """Return the branch shifts in lexicographic CRT-bit order."""

    from itertools import product

    prime_tuple = tuple(primes)
    return tuple(crt(bits, prime_tuple) for bits in product((0, 1), repeat=len(prime_tuple)))
