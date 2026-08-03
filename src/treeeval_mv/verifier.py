"""Primary exact verifiers."""

from __future__ import annotations

from dataclasses import dataclass

from .model import MVFamily, PolynomialMatchingFamily


@dataclass(frozen=True)
class VerificationReport:
    ok: bool
    checks: int
    errors: tuple[str, ...]


def dot_mod(left: tuple[int, ...], right: tuple[int, ...], modulus: int) -> int:
    if len(left) != len(right):
        raise ValueError("inner-product dimensions differ")
    return sum(a * b for a, b in zip(left, right, strict=True)) % modulus


def verify_mv_family(family: MVFamily, *, max_errors: int = 20) -> VerificationReport:
    errors: list[str] = []
    checks = 0
    if len(family.u) != len(family.v):
        errors.append(f"U has {len(family.u)} rows but V has {len(family.v)}")
        return VerificationReport(False, checks, tuple(errors))
    dimensions = {len(row) for row in family.u + family.v}
    if len(dimensions) > 1:
        errors.append(f"inconsistent vector dimensions: {sorted(dimensions)}")
        return VerificationReport(False, checks, tuple(errors))
    if family.diagonal_value in family.allowed_offdiagonal:
        errors.append("diagonal value is also allowed off diagonal")
    for i, left in enumerate(family.u):
        for j, right in enumerate(family.v):
            value = dot_mod(left, right, family.modulus)
            checks += 1
            if i == j and value != family.diagonal_value:
                errors.append(
                    f"diagonal ({i},{j}) is {value}, expected {family.diagonal_value}"
                )
            elif i != j and value not in family.allowed_offdiagonal:
                errors.append(
                    f"off diagonal ({i},{j}) is {value}, allowed={sorted(family.allowed_offdiagonal)}"
                )
            if len(errors) >= max_errors:
                return VerificationReport(False, checks, tuple(errors))
    return VerificationReport(not errors, checks, tuple(errors))


def evaluate_multilinear_polynomial(
    polynomial: dict[tuple[int, ...], int], point: tuple[int, ...], modulus: int
) -> int:
    total = 0
    for monomial, coefficient in polynomial.items():
        term = coefficient
        for variable in monomial:
            term *= point[variable]
        total += term
    return total % modulus


def verify_polynomial_matching_family(
    family: PolynomialMatchingFamily, *, max_errors: int = 20
) -> VerificationReport:
    errors: list[str] = []
    checks = 0
    if len(family.polynomials) != len(family.points):
        return VerificationReport(False, 0, ("polynomial and point counts differ",))
    basis_set = set(family.basis)
    for index, polynomial in enumerate(family.polynomials):
        outside = set(polynomial) - basis_set
        if outside:
            errors.append(f"polynomial {index} uses monomials outside the declared basis")
    for polynomial_index, polynomial in enumerate(family.polynomials):
        for point_index, point in enumerate(family.points):
            value = evaluate_multilinear_polynomial(polynomial, point, family.modulus)
            checks += 1
            if polynomial_index == point_index and value != 0:
                errors.append(
                    f"polynomial diagonal ({polynomial_index},{point_index}) is {value}, expected 0"
                )
            elif polynomial_index != point_index and value not in family.allowed_offdiagonal:
                errors.append(
                    f"polynomial off diagonal ({polynomial_index},{point_index}) is {value}"
                )
            if len(errors) >= max_errors:
                return VerificationReport(False, checks, tuple(errors))
    return VerificationReport(not errors, checks, tuple(errors))
