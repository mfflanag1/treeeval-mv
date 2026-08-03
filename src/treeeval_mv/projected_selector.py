"""Exact checks for a characteristic-local variant of the HPR selector.

HPR selects one nonzero monomial of a four-term polynomial at every CRT
prime.  At one distinguished prime, coefficient extraction can instead be
replaced by the polynomial's mixed finite difference.  This removes the
``(alpha, beta)`` state for that prime, but only for an output in the same
characteristic.  The routines here verify both the positive identity and the
cross-characteristic obstruction.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product

from .hpr_selector import _selector_monomial
from .model import MVFamily, canonical_hpr_offdiagonal_set


@dataclass(frozen=True)
class ProjectedSelectorReport:
    ok: bool
    brute_force_checks: int
    factor_checks: int
    errors: tuple[str, ...]


def _signed(value: int, parity: int, modulus: int) -> int:
    return (-value if parity & 1 else value) % modulus


def _projected_sum(
    *,
    primes: tuple[int, ...],
    distinguished: int,
    g1: tuple[int, ...],
    g2: tuple[int, ...],
    delta1: tuple[int, ...],
    delta2: tuple[int, ...],
) -> int:
    """Return the hybrid selector coefficient in the distinguished field."""

    output_prime = primes[distinguished]
    monomials: list[tuple[int, int] | None] = []
    normalizer = 1
    for index, prime in enumerate(primes):
        if index == distinguished:
            monomials.append(None)
            continue
        alpha, beta = _selector_monomial(prime, g1[index], g2[index])
        monomials.append((alpha, beta))
        normalizer = normalizer * (alpha % output_prime) % output_prime
    inverse_normalizer = pow(normalizer, -1, output_prime)

    total = 0
    bit_count = 2 * len(primes)
    for bits in product((0, 1), repeat=bit_count):
        left_bits = bits[: len(primes)]
        right_bits = bits[len(primes) :]
        include = True
        for index, prime in enumerate(primes):
            if index == distinguished:
                continue
            monomial = monomials[index]
            assert monomial is not None
            _, beta = monomial
            value = (
                (g1[index] + left_bits[index] * delta1[index])
                * (g2[index] + right_bits[index] * delta2[index])
            ) % prime
            if value != beta:
                include = False
                break
        if not include:
            continue
        raw = (
            (g1[distinguished] + left_bits[distinguished] * delta1[distinguished])
            * (g2[distinguished] + right_bits[distinguished] * delta2[distinguished])
        ) % output_prime
        parity = sum(left_bits) + sum(right_bits)
        total += _signed(raw * inverse_normalizer, parity, output_prime)
    return total % output_prime


def verify_projected_selector_identity(
    family: MVFamily, *, max_errors: int = 20
) -> ProjectedSelectorReport:
    """Exhaustively verify the hybrid selector for every CRT component.

    The check is universal over all possible base inner products and all
    canonical ``{0,1}`` delta patterns, rather than sampling family entries.
    """

    errors: list[str] = []
    brute_force_checks = 0
    factor_checks = 0
    primes = family.prime_factors
    if family.convention != "hpr" or family.diagonal_value != 1:
        return ProjectedSelectorReport(
            False, 0, 0, ("projected selector requires the HPR diagonal-one convention",)
        )
    if not primes or any(prime % 2 == 0 for prime in primes):
        return ProjectedSelectorReport(
            False, 0, 0, ("projected selector requires distinct odd prime factors",)
        )
    if family.allowed_offdiagonal != canonical_hpr_offdiagonal_set(primes):
        return ProjectedSelectorReport(
            False, 0, 0, ("family does not use HPR's canonical off-diagonal set",)
        )

    residue_ranges = tuple(range(prime) for prime in primes)
    delta_patterns = tuple(product((0, 1), repeat=len(primes)))
    for distinguished, output_prime in enumerate(primes):
        for g1 in product(*residue_ranges):
            for g2 in product(*residue_ranges):
                for delta1 in delta_patterns:
                    for delta2 in delta_patterns:
                        actual = _projected_sum(
                            primes=primes,
                            distinguished=distinguished,
                            g1=g1,
                            g2=g2,
                            delta1=delta1,
                            delta2=delta2,
                        )
                        expected = int(all(delta1) and all(delta2))
                        brute_force_checks += 1
                        if actual != expected:
                            errors.append(
                                f"prime {output_prime}: projected sum {actual}, "
                                f"expected {expected}; g=({g1},{g2}), "
                                f"delta=({delta1},{delta2})"
                            )
                            if len(errors) >= max_errors:
                                return ProjectedSelectorReport(
                                    False,
                                    brute_force_checks,
                                    factor_checks,
                                    tuple(errors),
                                )

                        # Independent factor check.  The distinguished factor
                        # is the mixed finite difference of xy; every other
                        # factor is HPR's signed coefficient indicator.
                        derivative = 0
                        for left_bit, right_bit in product((0, 1), repeat=2):
                            raw = (
                                (g1[distinguished] + left_bit * delta1[distinguished])
                                * (g2[distinguished] + right_bit * delta2[distinguished])
                            ) % output_prime
                            derivative += _signed(
                                raw, left_bit + right_bit, output_prime
                            )
                        derivative %= output_prime
                        factor_checks += 1
                        if derivative != delta1[distinguished] * delta2[distinguished]:
                            errors.append(
                                f"prime {output_prime}: mixed difference {derivative} "
                                "does not equal the product of deltas"
                            )
                            if len(errors) >= max_errors:
                                return ProjectedSelectorReport(
                                    False,
                                    brute_force_checks,
                                    factor_checks,
                                    tuple(errors),
                                )

    return ProjectedSelectorReport(
        not errors, brute_force_checks, factor_checks, tuple(errors)
    )


def has_cross_characteristic_fixed_weight(source_prime: int, output_prime: int) -> bool:
    """Exhaustively test for a forbidden fixed derivative weight on small primes.

    We ask whether some function ``phi: F_source -> F_output`` makes

    ``phi(xy)-phi((x+1)y)-phi(x(y+1))+phi((x+1)(y+1)) == 1``

    for every ``x,y``.  The general telescoping proof says no for distinct odd
    primes; exhaustive enumeration is retained as an independent finite check.
    """

    if source_prime < 2 or output_prime < 2:
        raise ValueError("prime arguments must be at least two")
    for values in product(range(output_prime), repeat=source_prime):
        works = True
        for x, y in product(range(source_prime), repeat=2):
            terms = (
                values[(x * y) % source_prime]
                - values[((x + 1) * y) % source_prime]
                - values[(x * (y + 1)) % source_prime]
                + values[((x + 1) * (y + 1)) % source_prime]
            ) % output_prime
            if terms != 1 % output_prime:
                works = False
                break
        if works:
            return True
    return False
