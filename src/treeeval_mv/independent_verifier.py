"""A deliberately separate CRT-based MV verifier.

This implementation does not import the primary verifier and computes scalar
products prime-by-prime using double-and-add multiplication.  It is intended as
an independent check before accepting any explicit family.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import prod

from .model import MVFamily, crt


@dataclass(frozen=True)
class IndependentReport:
    ok: bool
    checks: int
    errors: tuple[str, ...]


def _multiply_mod(left: int, right: int, modulus: int) -> int:
    result = 0
    addend = left % modulus
    multiplier = right % modulus
    while multiplier:
        if multiplier & 1:
            result = (result + addend) % modulus
        addend = (addend + addend) % modulus
        multiplier >>= 1
    return result


def _primewise_dot(left: tuple[int, ...], right: tuple[int, ...], primes: tuple[int, ...]) -> int:
    residues: list[int] = []
    for prime in primes:
        accumulator = 0
        for a, b in zip(left, right, strict=True):
            accumulator = (accumulator + _multiply_mod(a, b, prime)) % prime
        residues.append(accumulator)
    return crt(residues, primes)


def verify_mv_family_independently(
    family: MVFamily, *, max_errors: int = 20
) -> IndependentReport:
    errors: list[str] = []
    checks = 0
    primes = family.prime_factors
    if not primes or prod(primes) != family.modulus:
        return IndependentReport(
            False, 0, ("independent verifier requires a squarefree prime factorization",)
        )
    if len(family.u) != len(family.v):
        return IndependentReport(False, 0, ("U and V counts differ",))
    dimensions = [len(row) for row in family.u + family.v]
    if dimensions and any(value != dimensions[0] for value in dimensions):
        return IndependentReport(False, 0, ("vector dimensions differ",))
    for i in range(family.size):
        for j in range(family.size):
            value = _primewise_dot(family.u[i], family.v[j], primes)
            checks += 1
            valid = (
                value == family.diagonal_value
                if i == j
                else value in family.allowed_offdiagonal
            )
            if not valid:
                errors.append(f"entry ({i},{j}) has independently computed value {value}")
            if len(errors) >= max_errors:
                return IndependentReport(False, checks, tuple(errors))
    return IndependentReport(not errors, checks, tuple(errors))
