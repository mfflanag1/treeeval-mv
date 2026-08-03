"""Verifier for the polynomial selector identity used in HPR Lemmas 3.7--3.10."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product

from .model import MVFamily, canonical_hpr_offdiagonal_set


@dataclass(frozen=True)
class SelectorReport:
    ok: bool
    checks: int
    errors: tuple[str, ...]


def _dot_prime(left: tuple[int, ...], right: tuple[int, ...], prime: int) -> int:
    return sum((a % prime) * (b % prime) for a, b in zip(left, right, strict=True)) % prime


def _selector_monomial(prime: int, g1: int, g2: int) -> tuple[int, int]:
    coefficients: dict[int, int] = {}
    terms = (
        ((g1 * g2) % prime, 1),
        (((g1 + 1) * g2) % prime, -1),
        ((g1 * (g2 + 1)) % prime, -1),
        (((g1 + 1) * (g2 + 1)) % prime, 1),
    )
    for exponent, coefficient in terms:
        coefficients[exponent] = coefficients.get(exponent, 0) + coefficient
    nonzero = [(exponent, coefficient) for exponent, coefficient in coefficients.items() if coefficient]
    if not nonzero:
        raise AssertionError("HPR Lemma 3.7 polynomial vanished")
    beta, alpha = min(nonzero)
    if alpha not in {-2, -1, 1, 2}:
        raise AssertionError(f"unexpected selector coefficient {alpha}")
    return alpha, beta


def _signed_level_sum(
    *, prime: int, g1: int, g2: int, delta1: int, delta2: int, beta: int
) -> int:
    total = 0
    for bit1, bit2 in product((0, 1), repeat=2):
        value = ((g1 + bit1 * delta1) * (g2 + bit2 * delta2)) % prime
        if value == beta:
            total += -1 if (bit1 + bit2) & 1 else 1
    return total


def verify_hpr_selector_identity(
    family: MVFamily, *, max_errors: int = 20
) -> SelectorReport:
    """Exhaustively verify HPR's primewise selector and global isolation logic.

    This requires odd primes.  Modulus 6 is valid for the customary Grolmusz MV
    calibration but cannot calibrate HPR's selector because a coefficient ±2
    must be invertible modulo every prime factor.
    """

    errors: list[str] = []
    checks = 0
    primes = family.prime_factors
    if family.convention != "hpr" or family.diagonal_value != 1:
        return SelectorReport(False, 0, ("selector verifier requires the HPR diagonal-one convention",))
    if not primes or any(prime % 2 == 0 for prime in primes):
        return SelectorReport(False, 0, ("HPR selector requires distinct odd prime factors",))
    expected_offdiagonal = canonical_hpr_offdiagonal_set(primes)
    if family.allowed_offdiagonal != expected_offdiagonal:
        return SelectorReport(False, 0, ("family does not use HPR's canonical off-diagonal set",))

    # Lemma 3.9, including every possible g1,g2.  For an actual family entry,
    # delta is the corresponding inner product modulo p.
    for prime in primes:
        for g1, g2 in product(range(prime), repeat=2):
            alpha, beta = _selector_monomial(prime, g1, g2)
            if alpha % prime == 0:
                errors.append(f"selector alpha={alpha} is not a unit modulo {prime}")
            for a, b, r, s in product(range(family.size), repeat=4):
                delta1 = _dot_prime(family.u[a], family.v[r], prime)
                delta2 = _dot_prime(family.u[b], family.v[s], prime)
                total = _signed_level_sum(
                    prime=prime,
                    g1=g1,
                    g2=g2,
                    delta1=delta1,
                    delta2=delta2,
                    beta=beta,
                )
                checks += 1
                if r == a and s == b and total != alpha:
                    errors.append(
                        f"prime {prime}: target selector sum {total}, expected alpha {alpha}"
                    )
                elif delta1 * delta2 % prime == 0 and total != 0:
                    errors.append(
                        f"prime {prime}: zero-delta selector sum is {total}, expected zero"
                    )
                if len(errors) >= max_errors:
                    return SelectorReport(False, checks, tuple(errors))

    # Lemma 3.10's only family-dependent premise: every non-target pair has a
    # zero residue at some prime.  This is the exact isolation property used.
    for a, r in product(range(family.size), repeat=2):
        checks += 1
        residues = [_dot_prime(family.u[a], family.v[r], prime) for prime in primes]
        if a == r and any(value != 1 for value in residues):
            errors.append(f"diagonal ({a},{r}) does not have all-one CRT residues")
        elif a != r and all(value != 0 for value in residues):
            errors.append(f"off diagonal ({a},{r}) has no zero CRT residue")
        if len(errors) >= max_errors:
            return SelectorReport(False, checks, tuple(errors))
    return SelectorReport(not errors, checks, tuple(errors))
