"""Implicit coordinate-equality vectors and exact recursive-cost accounting.

This is a prototype changed vector type, not a TreeEval algorithm.  A label
word generates matching-vector coordinates on demand, and a single logical
overlay can be queried and reversed without materializing the update vector.
The accounting object makes the remaining nesting cost explicit.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil, log2, prod

from .model import crt


@dataclass(frozen=True)
class FactoredEqualityVectors:
    alphabet_size: int
    primes: tuple[int, ...]

    def __post_init__(self) -> None:
        if self.alphabet_size < 2:
            raise ValueError("alphabet_size must be at least two")
        if not self.primes or len(set(self.primes)) != len(self.primes):
            raise ValueError("primes must be nonempty and distinct")

    @property
    def modulus(self) -> int:
        return prod(self.primes)

    @property
    def capacity(self) -> int:
        return self.alphabet_size ** len(self.primes)

    def word(self, label: int) -> tuple[int, ...]:
        """Return the fixed-width base-``alphabet_size`` word for ``label``."""

        if not 0 <= label < self.capacity:
            raise ValueError("label is outside the family capacity")
        digits = [0] * len(self.primes)
        value = label
        for position in range(len(digits) - 1, -1, -1):
            value, digits[position] = divmod(value, self.alphabet_size)
        return tuple(digits)

    def coordinate(self, label: int, coordinate: int) -> int:
        if not 0 <= coordinate < self.alphabet_size:
            raise ValueError("coordinate is outside the vector dimension")
        word = self.word(label)
        return crt(
            [int(symbol == coordinate) for symbol in word],
            self.primes,
        )

    def materialize(self, label: int) -> tuple[int, ...]:
        return tuple(
            self.coordinate(label, coordinate)
            for coordinate in range(self.alphabet_size)
        )

    def gram(self, left_label: int, right_label: int) -> int:
        return sum(
            self.coordinate(left_label, coordinate)
            * self.coordinate(right_label, coordinate)
            for coordinate in range(self.alphabet_size)
        ) % self.modulus


@dataclass(frozen=True)
class Overlay:
    label: int
    scalar: int


def query_overlay_direct(
    vectors: FactoredEqualityVectors,
    base: tuple[int, ...],
    overlay: Overlay,
    query_label: int,
) -> int:
    """Query ``base + scalar*u_label`` by direct coordinate generation."""

    if len(base) != vectors.alphabet_size:
        raise ValueError("base vector has the wrong dimension")
    modulus = vectors.modulus
    return sum(
        (
            base[coordinate]
            + overlay.scalar * vectors.coordinate(overlay.label, coordinate)
        )
        * vectors.coordinate(query_label, coordinate)
        for coordinate in range(vectors.alphabet_size)
    ) % modulus


def query_overlay_primewise(
    vectors: FactoredEqualityVectors,
    base: tuple[int, ...],
    overlay: Overlay,
    query_label: int,
) -> int:
    """Independent query using component equality rather than CRT coordinates."""

    if len(base) != vectors.alphabet_size:
        raise ValueError("base vector has the wrong dimension")
    left_word = vectors.word(overlay.label)
    right_word = vectors.word(query_label)
    residues = []
    for component, prime in enumerate(vectors.primes):
        base_inner = base[right_word[component]] % prime
        delta = int(left_word[component] == right_word[component])
        residues.append((base_inner + overlay.scalar * delta) % prime)
    return crt(residues, vectors.primes)


def apply_overlay(
    vectors: FactoredEqualityVectors,
    base: tuple[int, ...],
    overlay: Overlay,
) -> tuple[int, ...]:
    if len(base) != vectors.alphabet_size:
        raise ValueError("base vector has the wrong dimension")
    return tuple(
        (base[coordinate] + overlay.scalar * vectors.coordinate(overlay.label, coordinate))
        % vectors.modulus
        for coordinate in range(vectors.alphabet_size)
    )


@dataclass(frozen=True)
class RecursiveTypeAccounting:
    height: int
    label_bits: int
    prime_count: int
    calls_per_level: int
    oracle_leaf_count: int
    branch_exponent: int
    suspended_overlay_bits: int
    full_overlay_stack_bits: int


def recursive_type_accounting(
    *, height: int, label_bits: int, prime_count: int, calls_per_level: int
) -> RecursiveTypeAccounting:
    """Account exactly for naive nesting of the implicit overlay type."""

    if min(height, label_bits, prime_count) < 0 or calls_per_level < 1:
        raise ValueError("accounting parameters are outside their valid ranges")
    suspended = label_bits + prime_count + ceil(log2(calls_per_level))
    branch_exponent = height * ceil(log2(calls_per_level))
    return RecursiveTypeAccounting(
        height=height,
        label_bits=label_bits,
        prime_count=prime_count,
        calls_per_level=calls_per_level,
        oracle_leaf_count=calls_per_level**height,
        branch_exponent=branch_exponent,
        suspended_overlay_bits=suspended,
        full_overlay_stack_bits=height * suspended,
    )
