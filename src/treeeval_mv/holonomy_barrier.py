"""Finite-state holonomy barrier for cross-characteristic selectors.

The mixed-difference equation with constant right-hand side cannot be a
single-valued potential on ``F_p^2`` with values in ``F_q`` when ``p`` and
``q`` are distinct.  A reversible stateful repair must remember the resulting
target-field holonomy.  This module records the exact state-count lower bound
and an independent finite permutation exhaustion for small parameters.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations
from math import gcd


@dataclass(frozen=True)
class HolonomyReport:
    source_modulus: int
    target_modulus: int
    cell_increment: int
    macro_increment: int
    minimum_states: int
    minimum_bits: int
    exhaustive_state_counts: tuple[tuple[int, bool], ...]
    construction_restores: bool


def minimum_reversible_holonomy_states(
    *, source_modulus: int, target_modulus: int, cell_increment: int = 1
) -> int:
    """Return the minimum cycle length that can absorb the target holonomy.

    Summing a constant cell increment around one source cycle produces
    ``source_modulus * cell_increment`` in the target group.  If a reversible
    state returns after ``L`` source cycles, exact restoration requires the
    target modulus to divide ``L`` times that macro increment.
    """

    if source_modulus < 2 or target_modulus < 2:
        raise ValueError("moduli must be at least two")
    macro_increment = source_modulus * cell_increment
    return target_modulus // gcd(macro_increment, target_modulus)


def minimum_state_bits(state_count: int) -> int:
    if state_count < 1:
        raise ValueError("state_count must be positive")
    return (state_count - 1).bit_length()


def _permutation_cycles(permutation: tuple[int, ...]) -> tuple[int, ...]:
    seen: set[int] = set()
    lengths: list[int] = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            current = permutation[current]
            length += 1
        lengths.append(length)
    return tuple(lengths)


def has_reversible_holonomy_permutation(
    *,
    state_count: int,
    source_modulus: int,
    target_modulus: int,
    cell_increment: int = 1,
) -> bool:
    """Exhaust permutations whose every cycle restores arbitrary state."""

    if state_count < 1:
        return False
    macro_increment = (source_modulus * cell_increment) % target_modulus
    for permutation in permutations(range(state_count)):
        if all(
            length * macro_increment % target_modulus == 0
            for length in _permutation_cycles(permutation)
        ):
            return True
    return False


def verify_holonomy_bound(
    *, source_modulus: int, target_modulus: int, cell_increment: int = 1
) -> HolonomyReport:
    """Cross-check the number-theoretic bound by finite permutation search."""

    minimum = minimum_reversible_holonomy_states(
        source_modulus=source_modulus,
        target_modulus=target_modulus,
        cell_increment=cell_increment,
    )
    exhaustive = tuple(
        (
            states,
            has_reversible_holonomy_permutation(
                state_count=states,
                source_modulus=source_modulus,
                target_modulus=target_modulus,
                cell_increment=cell_increment,
            ),
        )
        for states in range(1, minimum + 1)
    )
    macro_increment = (source_modulus * cell_increment) % target_modulus
    construction_restores = minimum * macro_increment % target_modulus == 0
    return HolonomyReport(
        source_modulus=source_modulus,
        target_modulus=target_modulus,
        cell_increment=cell_increment,
        macro_increment=macro_increment,
        minimum_states=minimum,
        minimum_bits=minimum_state_bits(minimum),
        exhaustive_state_counts=exhaustive,
        construction_restores=construction_restores,
    )
