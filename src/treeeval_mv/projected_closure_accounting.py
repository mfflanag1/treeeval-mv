"""Exact bit and recurrence ledger for composing projected selectors."""

from __future__ import annotations

from dataclasses import dataclass
from math import prod


def ceiling_log2(value: int) -> int:
    if value < 1:
        raise ValueError("value must be positive")
    return (value - 1).bit_length()


@dataclass(frozen=True)
class ProjectedClosureAccounting:
    primes: tuple[int, ...]
    modulus: int
    depth: int
    label_bits: int
    hpr_modulus_bits_per_level: int
    projected_complement_bits: tuple[int, ...]
    sequential_worst_bits_per_level: int
    hpr_stack_bits: int
    sequential_projected_stack_bits: int
    candidate_outer_oracle_log2_lower_bound: int


def projected_closure_accounting(
    *, primes: tuple[int, ...], depth: int, label_bits: int
) -> ProjectedClosureAccounting:
    """Account for the obvious full-update compositions.

    Sequentially running every characteristic-local output routine cannot use
    the best component's state bound: a recursive execution path eventually
    enters the routine with the largest complementary modulus.  Reordering by
    candidate avoids that stack but repeats a child-oracle shift for at least
    every candidate pair, giving ``N^(2h)`` calls for ``N=2^ell``.
    """

    if not primes or len(set(primes)) != len(primes):
        raise ValueError("primes must be a nonempty tuple of distinct values")
    if depth < 0 or label_bits < 0:
        raise ValueError("depth and label_bits must be nonnegative")
    modulus = prod(primes)
    complement_bits = tuple(ceiling_log2(modulus // prime) for prime in primes)
    worst = max(complement_bits)
    return ProjectedClosureAccounting(
        primes=primes,
        modulus=modulus,
        depth=depth,
        label_bits=label_bits,
        hpr_modulus_bits_per_level=ceiling_log2(modulus),
        projected_complement_bits=complement_bits,
        sequential_worst_bits_per_level=worst,
        hpr_stack_bits=depth * ceiling_log2(modulus),
        sequential_projected_stack_bits=depth * worst,
        candidate_outer_oracle_log2_lower_bound=2 * label_bits * depth,
    )
