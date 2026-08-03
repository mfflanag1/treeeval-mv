"""Joint dense packing for arbitrary-initialized Z_m vector registers.

HPR Remark 2.3 validates each coordinate separately.  If all D coordinates
are packed as one integer in [0,m**D), only one quotient bit is needed to map
an arbitrary b=ceil(log2(m**D))-bit catalyst into the valid range.  The helpers
below implement the exact state map used by the accompanying packing lemma.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence


def packed_bit_length(modulus: int, coordinate_count: int) -> int:
    if modulus < 2 or coordinate_count < 1:
        raise ValueError("modulus >= 2 and coordinate_count >= 1 are required")
    state_count = modulus**coordinate_count
    return (state_count - 1).bit_length()


def pack_digits(digits: Sequence[int], modulus: int) -> int:
    value = 0
    place = 1
    for digit in digits:
        if not 0 <= digit < modulus:
            raise ValueError(f"digit {digit} is outside [0,{modulus})")
        value += digit * place
        place *= modulus
    return value


def unpack_digits(value: int, modulus: int, coordinate_count: int) -> tuple[int, ...]:
    limit = modulus**coordinate_count
    if not 0 <= value < limit:
        raise ValueError(f"packed value must be in [0,{limit})")
    digits: list[int] = []
    remaining = value
    for _ in range(coordinate_count):
        remaining, digit = divmod(remaining, modulus)
        digits.append(digit)
    return tuple(digits)


@dataclass(frozen=True)
class NormalizedCatalyst:
    packed_vector: int
    quotient_bit: int
    bit_length: int
    state_count: int


def normalize_arbitrary_catalyst(
    bit_string_value: int, *, modulus: int, coordinate_count: int
) -> NormalizedCatalyst:
    """Map an arbitrary b-bit value to a valid packed vector plus one bit."""

    bits = packed_bit_length(modulus, coordinate_count)
    if not 0 <= bit_string_value < 1 << bits:
        raise ValueError(f"input is not a {bits}-bit value")
    state_count = modulus**coordinate_count
    quotient, remainder = divmod(bit_string_value, state_count)
    # Minimal b guarantees 2**b < 2*state_count, unless state_count is a
    # power of two (where the quotient is always zero).
    if quotient not in (0, 1):
        raise AssertionError("minimal dense packing unexpectedly needs more than one quotient bit")
    return NormalizedCatalyst(remainder, quotient, bits, state_count)


def restore_arbitrary_catalyst(normalized: NormalizedCatalyst) -> int:
    value = normalized.packed_vector + normalized.quotient_bit * normalized.state_count
    if value >= 1 << normalized.bit_length:
        raise AssertionError("restored value does not fit the original tape")
    return value


def update_coordinate(
    packed_vector: int,
    *,
    coordinate: int,
    delta: int,
    modulus: int,
    coordinate_count: int,
) -> int:
    """Add delta to one base-m digit without unpacking the stored state."""

    if not 0 <= coordinate < coordinate_count:
        raise ValueError("coordinate is out of range")
    limit = modulus**coordinate_count
    if not 0 <= packed_vector < limit:
        raise ValueError("packed vector is invalid")
    place = modulus**coordinate
    old_digit = (packed_vector // place) % modulus
    new_digit = (old_digit + delta) % modulus
    updated = packed_vector + (new_digit - old_digit) * place
    if not 0 <= updated < limit:
        raise AssertionError("coordinate update left the packed-vector range")
    return updated


def add_vector(
    packed_vector: int, vector: Iterable[int], *, modulus: int, coordinate_count: int
) -> int:
    result = packed_vector
    materialized = tuple(vector)
    if len(materialized) != coordinate_count:
        raise ValueError("vector length differs from packed coordinate count")
    for coordinate, delta in enumerate(materialized):
        result = update_coordinate(
            result,
            coordinate=coordinate,
            delta=delta,
            modulus=modulus,
            coordinate_count=coordinate_count,
        )
    return result


def inner_product(
    packed_vector: int, vector: Iterable[int], *, modulus: int, coordinate_count: int
) -> int:
    materialized = tuple(vector)
    if len(materialized) != coordinate_count:
        raise ValueError("vector length differs from packed coordinate count")
    total = 0
    place = 1
    for coefficient in materialized:
        digit = (packed_vector // place) % modulus
        total = (total + digit * coefficient) % modulus
        place *= modulus
    return total


@dataclass(frozen=True)
class PackingExhaustionReport:
    ok: bool
    initial_states: int
    update_checks: int
    errors: tuple[str, ...]


def exhaust_packing_state_map(
    *, modulus: int, coordinate_count: int, max_errors: int = 20
) -> PackingExhaustionReport:
    """Exhaust normalization, every single-coordinate update, and restoration."""

    errors: list[str] = []
    bits = packed_bit_length(modulus, coordinate_count)
    update_checks = 0
    for initial in range(1 << bits):
        normalized = normalize_arbitrary_catalyst(
            initial, modulus=modulus, coordinate_count=coordinate_count
        )
        if restore_arbitrary_catalyst(normalized) != initial:
            errors.append(f"initial state {initial} failed immediate restoration")
        original_vector = normalized.packed_vector
        for coordinate in range(coordinate_count):
            for delta in range(modulus):
                updated = update_coordinate(
                    original_vector,
                    coordinate=coordinate,
                    delta=delta,
                    modulus=modulus,
                    coordinate_count=coordinate_count,
                )
                restored_vector = update_coordinate(
                    updated,
                    coordinate=coordinate,
                    delta=-delta,
                    modulus=modulus,
                    coordinate_count=coordinate_count,
                )
                restored = restore_arbitrary_catalyst(
                    NormalizedCatalyst(
                        restored_vector,
                        normalized.quotient_bit,
                        normalized.bit_length,
                        normalized.state_count,
                    )
                )
                update_checks += 1
                if restored != initial:
                    errors.append(
                        f"state {initial}, coordinate {coordinate}, delta {delta} failed restoration"
                    )
                if len(errors) >= max_errors:
                    return PackingExhaustionReport(
                        False, 1 << bits, update_checks, tuple(errors)
                    )
    return PackingExhaustionReport(not errors, 1 << bits, update_checks, tuple(errors))
