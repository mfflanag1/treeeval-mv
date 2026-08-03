"""Data structures and convention conversions for matching-vector families.

There are two conventions in the literature.  The customary DGY convention has
diagonal inner product zero and nonzero off-diagonal products.  HPR uses the
affinely flipped convention: diagonal one, every prime residue is in {0, 1},
and an off-diagonal product is never one.  Keeping these conventions explicit
prevents a very easy but consequential verifier error.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import product
from math import prod
from typing import Any, Iterable, Mapping, Sequence

Vector = tuple[int, ...]
Monomial = tuple[int, ...]
Polynomial = dict[Monomial, int]


def _normalized_vectors(rows: Iterable[Sequence[int]], modulus: int) -> tuple[Vector, ...]:
    return tuple(tuple(int(x) % modulus for x in row) for row in rows)


def crt(residues: Sequence[int], moduli: Sequence[int]) -> int:
    """Return the Chinese-remainder representative in [0, product(moduli))."""

    if len(residues) != len(moduli) or not moduli:
        raise ValueError("residues and moduli must have equal nonzero length")
    modulus = prod(moduli)
    answer = 0
    for residue, prime in zip(residues, moduli, strict=True):
        partial = modulus // prime
        answer += (residue % prime) * partial * pow(partial, -1, prime)
    return answer % modulus


def canonical_cube(primes: Sequence[int]) -> frozenset[int]:
    """All CRT values having residue zero or one at every supplied prime."""

    if len(set(primes)) != len(primes) or any(p < 2 for p in primes):
        raise ValueError("primes must be distinct integers at least two")
    return frozenset(crt(bits, primes) for bits in product((0, 1), repeat=len(primes)))


def canonical_customary_set(primes: Sequence[int]) -> frozenset[int]:
    """DGY's canonical off-diagonal set: the {0,1} CRT cube minus zero."""

    return canonical_cube(primes) - {0}


def canonical_hpr_offdiagonal_set(primes: Sequence[int]) -> frozenset[int]:
    """HPR's possible off-diagonal products: the CRT cube minus one."""

    return canonical_cube(primes) - {1}


@dataclass(frozen=True)
class MVFamily:
    modulus: int
    u: tuple[Vector, ...]
    v: tuple[Vector, ...]
    allowed_offdiagonal: frozenset[int]
    diagonal_value: int = 0
    prime_factors: tuple[int, ...] = ()
    convention: str = "customary"
    metadata: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def build(
        cls,
        *,
        modulus: int,
        u: Iterable[Sequence[int]],
        v: Iterable[Sequence[int]],
        allowed_offdiagonal: Iterable[int],
        diagonal_value: int = 0,
        prime_factors: Sequence[int] = (),
        convention: str = "customary",
        metadata: Mapping[str, Any] | None = None,
    ) -> "MVFamily":
        if modulus < 2:
            raise ValueError("modulus must be at least two")
        return cls(
            modulus=modulus,
            u=_normalized_vectors(u, modulus),
            v=_normalized_vectors(v, modulus),
            allowed_offdiagonal=frozenset(x % modulus for x in allowed_offdiagonal),
            diagonal_value=diagonal_value % modulus,
            prime_factors=tuple(prime_factors),
            convention=convention,
            metadata={} if metadata is None else dict(metadata),
        )

    @property
    def size(self) -> int:
        return len(self.u)

    @property
    def dimension(self) -> int:
        return len(self.u[0]) if self.u else 0


@dataclass(frozen=True)
class PolynomialMatchingFamily:
    modulus: int
    polynomials: tuple[Polynomial, ...]
    points: tuple[Vector, ...]
    basis: tuple[Monomial, ...]
    allowed_offdiagonal: frozenset[int]
    prime_factors: tuple[int, ...] = ()
    metadata: Mapping[str, Any] = field(default_factory=dict)

    @property
    def size(self) -> int:
        return len(self.polynomials)

    @property
    def dimension(self) -> int:
        return len(self.basis)


def flip_customary_to_hpr(family: MVFamily) -> MVFamily:
    """Apply HPR Remark 2.11: u'=(u,1), v'=(-v,1)."""

    if family.diagonal_value != 0:
        raise ValueError("the source family must use diagonal zero")
    modulus = family.modulus
    u = [row + (1,) for row in family.u]
    v = [tuple((-x) % modulus for x in row) + (1,) for row in family.v]
    return MVFamily.build(
        modulus=modulus,
        u=u,
        v=v,
        allowed_offdiagonal={(1 - s) % modulus for s in family.allowed_offdiagonal},
        diagonal_value=1,
        prime_factors=family.prime_factors,
        convention="hpr",
        metadata={**family.metadata, "affine_flip": "HPR Remark 2.11"},
    )
