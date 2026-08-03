"""Exact closure harness for stateless four-call selector kernels.

The HPR selector makes four shifted child calls.  The most direct way to
remove its primewise coefficient state is to apply one fixed scalar lookup
``phi`` to the product inner product returned at each shift.  Such a kernel is
recursively closed: if its one-level coefficient is the Kronecker delta, then
two levels compose by multiplying the two coefficients and use sixteen leaf
calls without retaining a modulus-sized selector descriptor.

This module checks that identity exactly and records a three-equation
obstruction for every modulus with a nontrivial CRT idempotent.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from math import prod
from typing import Sequence

from .model import crt


@dataclass(frozen=True)
class FixedWeightClosureReport:
    ok: bool
    modulus: int
    one_level_checks: int
    two_level_checks: int
    leaf_calls_per_two_levels: int
    errors: tuple[str, ...]


@dataclass(frozen=True)
class FixedWeightBarrierReport:
    impossible: bool
    modulus: int
    idempotent: int
    complement: int
    equations: tuple[str, ...]
    contradiction: str


def binary_crt_deltas(primes: Sequence[int]) -> tuple[int, ...]:
    """Return the residues whose CRT coordinates all lie in ``{0,1}``."""

    prime_tuple = tuple(primes)
    if not prime_tuple or prod(prime_tuple) < 2:
        raise ValueError("a nontrivial prime factorization is required")
    return tuple(
        crt(residues, prime_tuple)
        for residues in product((0, 1), repeat=len(prime_tuple))
    )


def mixed_difference(
    phi: Sequence[int], *, g1: int, g2: int, delta1: int, delta2: int
) -> int:
    """Evaluate the four-call scalar coefficient using exact ring arithmetic."""

    modulus = len(phi)
    if modulus < 2:
        raise ValueError("phi must contain one value for every residue")
    terms = (
        phi[(g1 * g2) % modulus],
        -phi[((g1 + delta1) * g2) % modulus],
        -phi[(g1 * (g2 + delta2)) % modulus],
        phi[((g1 + delta1) * (g2 + delta2)) % modulus],
    )
    return sum(terms) % modulus


def verify_fixed_weight_two_level_closure(
    phi: Sequence[int], *, primes: Sequence[int], exhaustive_two_level: bool = True
) -> FixedWeightClosureReport:
    """Verify one-level isolation and its symbolic two-level composition.

    A two-level leaf coefficient factors as the product of the two one-level
    coefficients.  Exhausting ordered pairs of one-level contexts therefore
    checks the exact nested oracle trace without enumerating truth tables.
    """

    prime_tuple = tuple(primes)
    modulus = prod(prime_tuple)
    if len(phi) != modulus:
        raise ValueError("phi length must equal the product of the primes")
    deltas = binary_crt_deltas(prime_tuple)
    contexts: list[tuple[int, int]] = []
    errors: list[str] = []
    one_level_checks = 0
    for g1, g2, delta1, delta2 in product(
        range(modulus), range(modulus), deltas, deltas
    ):
        actual = mixed_difference(
            phi, g1=g1, g2=g2, delta1=delta1, delta2=delta2
        )
        expected = int(delta1 == 1 and delta2 == 1)
        contexts.append((actual, expected))
        one_level_checks += 1
        if actual != expected and len(errors) < 20:
            errors.append(
                f"one level: g=({g1},{g2}), delta=({delta1},{delta2}), "
                f"actual={actual}, expected={expected}"
            )

    two_level_checks = 0
    if exhaustive_two_level:
        for outer_actual, outer_expected in contexts:
            for inner_actual, inner_expected in contexts:
                two_level_checks += 1
                actual = outer_actual * inner_actual % modulus
                expected = outer_expected * inner_expected
                if actual != expected and len(errors) < 20:
                    errors.append(
                        f"two levels: coefficient {actual}, expected {expected}"
                    )
    return FixedWeightClosureReport(
        not errors,
        modulus,
        one_level_checks,
        two_level_checks,
        16,
        tuple(errors),
    )


def fixed_weight_crt_barrier(primes: Sequence[int]) -> FixedWeightBarrierReport:
    """Give a three-instance contradiction for a composite CRT modulus.

    Let ``e`` be a nontrivial CRT idempotent and ``f=1-e``.  At ``g1=0`` and
    ``delta1=1``, the four-call equation becomes
    ``phi(g2+delta2)-phi(g2)``.  The non-target instances
    ``(g2,delta2)=(1,e)`` and ``(1+e,f)`` force the endpoints of a path from
    1 to 2 to have equal weights, while the target instance ``(1,1)`` forces
    their difference to be one.
    """

    prime_tuple = tuple(primes)
    if len(prime_tuple) < 2:
        raise ValueError("the CRT barrier requires at least two factors")
    modulus = prod(prime_tuple)
    e = crt((1,) + (0,) * (len(prime_tuple) - 1), prime_tuple)
    f = (1 - e) % modulus
    if e in (0, 1) or f in (0, 1):
        raise AssertionError("the selected CRT idempotents are not nontrivial")
    return FixedWeightBarrierReport(
        True,
        modulus,
        e,
        f,
        (
            f"phi({(1 + e) % modulus}) - phi(1) = 0",
            "phi(2) - phi(1) = 1",
            f"phi(2) - phi({(1 + e) % modulus}) = 0",
        ),
        "adding the two non-target equations gives phi(2)-phi(1)=0, "
        "contradicting the target equation",
    )
