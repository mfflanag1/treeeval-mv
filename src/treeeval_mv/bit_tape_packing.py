"""Bit-level transition model for the joint packed-catalyst lemma.

Long constants are supplied by bit generators.  In the formal simulation those
bits come from logspace-uniform powering/multiplication circuits; this module
checks that the only destructive tape operations needed are single-pass binary
addition and subtraction with one carry/borrow bit.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from .packed_catalyst import packed_bit_length

BitGenerator = Callable[[int], int]


@dataclass
class BitTape:
    """Fixed-width little-endian catalytic bit tape."""

    bits: list[int]

    @classmethod
    def from_int(cls, value: int, width: int) -> "BitTape":
        if width < 1 or not 0 <= value < 1 << width:
            raise ValueError("value does not fit the requested positive width")
        return cls([(value >> index) & 1 for index in range(width)])

    @property
    def width(self) -> int:
        return len(self.bits)

    def to_int(self) -> int:
        return sum(bit << index for index, bit in enumerate(self.bits))

    def snapshot(self) -> tuple[int, ...]:
        return tuple(self.bits)


def integer_bit_generator(value: int) -> BitGenerator:
    if value < 0:
        raise ValueError("bit generators represent nonnegative integers")
    return lambda index: (value >> index) & 1


def compare_with_generated_constant(
    tape: BitTape, constant_bit: BitGenerator, *, constant_extra_bit: int = 0
) -> int:
    """Return -1, 0, or 1 according as tape is below, equal to, or above C.

    `constant_extra_bit` is the bit of C immediately above the tape width.  It
    handles the power-of-two case C=2**width without allocating another cell.
    """

    if constant_extra_bit not in (0, 1):
        raise ValueError("constant_extra_bit must be Boolean")
    if constant_extra_bit:
        return -1
    for index in range(tape.width - 1, -1, -1):
        constant = constant_bit(index)
        if tape.bits[index] < constant:
            return -1
        if tape.bits[index] > constant:
            return 1
    return 0


def subtract_generated_constant_in_place(tape: BitTape, constant_bit: BitGenerator) -> None:
    """Compute tape <- tape-C, assuming the result is nonnegative."""

    borrow = 0
    for index in range(tape.width):
        difference = tape.bits[index] - constant_bit(index) - borrow
        if difference < 0:
            difference += 2
            borrow = 1
        else:
            borrow = 0
        tape.bits[index] = difference
    if borrow:
        raise ValueError("subtraction underflow")


def add_generated_constant_in_place(tape: BitTape, constant_bit: BitGenerator) -> None:
    """Compute tape <- tape+C, assuming the result fits the fixed width."""

    carry = 0
    for index in range(tape.width):
        total = tape.bits[index] + constant_bit(index) + carry
        tape.bits[index] = total & 1
        carry = total >> 1
    if carry:
        raise ValueError("addition overflow")


@dataclass(frozen=True)
class BitTapeNormalization:
    quotient_bit: int
    state_count: int
    original_snapshot: tuple[int, ...]


def normalize_bit_tape_in_place(
    tape: BitTape, *, modulus: int, coordinate_count: int
) -> BitTapeNormalization:
    expected_width = packed_bit_length(modulus, coordinate_count)
    if tape.width != expected_width:
        raise ValueError(f"expected a {expected_width}-bit tape")
    original = tape.snapshot()
    state_count = modulus**coordinate_count
    generator = integer_bit_generator(state_count)
    extra = (state_count >> tape.width) & 1
    quotient = int(
        compare_with_generated_constant(
            tape, generator, constant_extra_bit=extra
        )
        >= 0
    )
    if quotient:
        subtract_generated_constant_in_place(tape, generator)
    if tape.to_int() >= state_count:
        raise AssertionError("normalization did not produce a valid packed vector")
    return BitTapeNormalization(quotient, state_count, original)


def restore_bit_tape_in_place(tape: BitTape, normalization: BitTapeNormalization) -> None:
    if normalization.quotient_bit:
        add_generated_constant_in_place(
            tape, integer_bit_generator(normalization.state_count)
        )
    if tape.snapshot() != normalization.original_snapshot:
        raise AssertionError("bit-level catalyst restoration failed")


def update_base_m_digit_in_place(
    tape: BitTape,
    *,
    coordinate: int,
    delta: int,
    modulus: int,
    coordinate_count: int,
) -> None:
    """Bit-level destructive update, with digit extraction as the TC0 oracle.

    The Python model evaluates the digit directly.  The proof replaces this
    line by logspace evaluation of the uniform division circuit, storing only
    the O(log m)-bit digit before the destructive ripple-carry pass.
    """

    if not 0 <= coordinate < coordinate_count:
        raise ValueError("coordinate is out of range")
    value = tape.to_int()
    state_count = modulus**coordinate_count
    if value >= state_count:
        raise ValueError("tape is not a valid packed base-m vector")
    place = modulus**coordinate
    old_digit = (value // place) % modulus
    new_digit = (old_digit + delta) % modulus
    signed_change = (new_digit - old_digit) * place
    if signed_change >= 0:
        add_generated_constant_in_place(tape, integer_bit_generator(signed_change))
    else:
        subtract_generated_constant_in_place(
            tape, integer_bit_generator(-signed_change)
        )
    if tape.to_int() >= state_count:
        raise AssertionError("digit update left the packed-vector range")


@dataclass(frozen=True)
class BitTapeExhaustionReport:
    ok: bool
    initial_states: int
    transitions: int
    errors: tuple[str, ...]


def exhaust_bit_tape_transitions(
    *, modulus: int, coordinate_count: int, max_errors: int = 20
) -> BitTapeExhaustionReport:
    width = packed_bit_length(modulus, coordinate_count)
    transitions = 0
    errors: list[str] = []
    for initial in range(1 << width):
        for coordinate in range(coordinate_count):
            for delta in range(modulus):
                tape = BitTape.from_int(initial, width)
                normalization = normalize_bit_tape_in_place(
                    tape, modulus=modulus, coordinate_count=coordinate_count
                )
                update_base_m_digit_in_place(
                    tape,
                    coordinate=coordinate,
                    delta=delta,
                    modulus=modulus,
                    coordinate_count=coordinate_count,
                )
                update_base_m_digit_in_place(
                    tape,
                    coordinate=coordinate,
                    delta=-delta,
                    modulus=modulus,
                    coordinate_count=coordinate_count,
                )
                try:
                    restore_bit_tape_in_place(tape, normalization)
                except AssertionError as error:
                    errors.append(
                        f"initial={initial}, coordinate={coordinate}, delta={delta}: {error}"
                    )
                transitions += 1
                if len(errors) >= max_errors:
                    return BitTapeExhaustionReport(
                        False, 1 << width, transitions, tuple(errors)
                    )
    return BitTapeExhaustionReport(not errors, 1 << width, transitions, tuple(errors))
