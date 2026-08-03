"""Exact CRT Gram-system characterization and bounded SMT search.

An HPR-convention matching-vector family over a squarefree modulus is
equivalent to a collection of binary Gram matrices, one over each prime
field.  Every component matrix has rank at most the vector dimension, its
diagonal is one, and the entrywise product of the component matrices is the
identity matrix.  Conversely, factorizations of such matrices can be combined
coordinatewise by CRT to recover an HPR family.

The SMT encoding below searches the factors directly.  A separate verifier
recomputes the Gram matrices and their ranks from the resulting explicit
family; it does not trust the solver's claimed factorization.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from itertools import product as cartesian_product
from math import prod
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time
from functools import reduce
from typing import Callable, Sequence

from .model import MVFamily, canonical_hpr_offdiagonal_set, crt


Matrix = tuple[tuple[int, ...], ...]


@dataclass(frozen=True)
class GramCharacterizationReport:
    ok: bool
    component_grams: tuple[Matrix, ...]
    component_ranks: tuple[int, ...]
    checks: int
    errors: tuple[str, ...]


@dataclass(frozen=True)
class GramSMTResult:
    status: str
    size: int
    dimension: int
    primes: tuple[int, ...]
    encoding_sha256: str
    solver: str
    wall_seconds: float
    variables: int
    constraints: int
    solver_seed: int
    family: MVFamily | None
    verification: GramCharacterizationReport | None


@dataclass(frozen=True)
class ExtensionResult:
    extended: bool
    family: MVFamily | None
    component_candidates: tuple[int, ...]
    combinations_checked: int


def equality_product_hpr_family(
    *, primes: Sequence[int], alphabet_size: int
) -> MVFamily:
    """Construct the characteristic-free ``N=d^t`` product calibration.

    A label is a ``t``-tuple over an alphabet of size ``d``.  Component ``i``
    uses the equality Gram matrix on the label's ``i``th symbol.  The
    Hadamard product is therefore the identity.
    """

    prime_tuple = tuple(primes)
    if len(prime_tuple) < 2 or len(set(prime_tuple)) != len(prime_tuple):
        raise ValueError("at least two distinct primes are required")
    if alphabet_size < 1:
        raise ValueError("alphabet_size must be positive")
    labels = tuple(cartesian_product(range(alphabet_size), repeat=len(prime_tuple)))
    rows: list[tuple[int, ...]] = []
    for label in labels:
        rows.append(
            tuple(
                crt(
                    [int(label[component] == coordinate) for component in range(len(prime_tuple))],
                    prime_tuple,
                )
                for coordinate in range(alphabet_size)
            )
        )
    modulus = prod(prime_tuple)
    return MVFamily.build(
        modulus=modulus,
        u=rows,
        v=rows,
        allowed_offdiagonal=canonical_hpr_offdiagonal_set(prime_tuple),
        diagonal_value=1,
        prime_factors=prime_tuple,
        convention="hpr",
        metadata={
            "construction": "coordinate-equality product Gram system",
            "label_alphabet_size": alphabet_size,
        },
    )


def determinant_lift_is_safe(*, prime: int, rank_bound: int) -> bool:
    """Test the exact 0/1 determinant condition for lifting rank to Q.

    The maximal determinant inequality for an ``n`` by ``n`` zero-one matrix
    is ``|det| <= (n+1)^((n+1)/2) / 2^n``.  With ``n=rank_bound+1``, a prime
    strictly above that quantity cannot divide a nonzero minor.  Squaring the
    comparison keeps this predicate exact and avoids floating-point arithmetic.
    """

    if prime < 2 or rank_bound < 0:
        raise ValueError("prime must be at least two and rank_bound nonnegative")
    scaled_prime = prime * 2 ** (rank_bound + 1)
    return scaled_prime * scaled_prime > (rank_bound + 2) ** (rank_bound + 2)


def characteristic_free_size_bound(*, primes: Sequence[int], rank_bound: int) -> int | None:
    """Return ``d^t`` when every component rank safely lifts to Q."""

    prime_tuple = tuple(primes)
    if all(
        determinant_lift_is_safe(prime=prime, rank_bound=rank_bound)
        for prime in prime_tuple
    ):
        return rank_bound ** len(prime_tuple)
    return None


def matrix_rank_mod(matrix: Sequence[Sequence[int]], prime: int) -> int:
    """Compute exact row rank over ``F_prime`` by modular elimination."""

    if prime < 2:
        raise ValueError("prime must be at least two")
    if not matrix:
        return 0
    width = len(matrix[0])
    if any(len(row) != width for row in matrix):
        raise ValueError("matrix is ragged")
    work = [[value % prime for value in row] for row in matrix]
    pivot_row = 0
    for column in range(width):
        pivot = next(
            (row for row in range(pivot_row, len(work)) if work[row][column]), None
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        inverse = pow(work[pivot_row][column], -1, prime)
        work[pivot_row] = [(value * inverse) % prime for value in work[pivot_row]]
        for row in range(len(work)):
            if row == pivot_row or not work[row][column]:
                continue
            multiplier = work[row][column]
            work[row] = [
                (left - multiplier * right) % prime
                for left, right in zip(work[row], work[pivot_row], strict=True)
            ]
        pivot_row += 1
        if pivot_row == len(work):
            break
    return pivot_row


def verify_gram_characterization(family: MVFamily) -> GramCharacterizationReport:
    """Independently verify the binary-rank/Hadamard characterization."""

    errors: list[str] = []
    checks = 0
    primes = family.prime_factors
    if not primes or prod(primes) != family.modulus:
        return GramCharacterizationReport(
            False, (), (), 0, ("a squarefree prime factorization is required",)
        )
    if family.diagonal_value != 1 or family.convention != "hpr":
        errors.append("family must use the HPR diagonal-one convention")
    if len(family.u) != len(family.v):
        errors.append("U and V list lengths differ")
        return GramCharacterizationReport(False, (), (), checks, tuple(errors))
    size = family.size
    dimension = family.dimension
    if any(len(row) != dimension for row in family.u + family.v):
        errors.append("vector dimensions differ")
        return GramCharacterizationReport(False, (), (), checks, tuple(errors))

    grams: list[Matrix] = []
    ranks: list[int] = []
    for prime in primes:
        matrix = tuple(
            tuple(
                sum(
                    left * right
                    for left, right in zip(family.u[i], family.v[j], strict=True)
                )
                % prime
                for j in range(size)
            )
            for i in range(size)
        )
        grams.append(matrix)
        rank = matrix_rank_mod(matrix, prime)
        ranks.append(rank)
        if rank > dimension:
            errors.append(f"component modulo {prime} has rank {rank} > {dimension}")
        for i in range(size):
            for j in range(size):
                checks += 1
                if matrix[i][j] not in (0, 1):
                    errors.append(
                        f"component modulo {prime} entry ({i},{j}) is {matrix[i][j]}, not binary"
                    )
                if i == j and matrix[i][j] != 1:
                    errors.append(f"component modulo {prime} diagonal ({i},{i}) is not one")

    for i in range(size):
        for j in range(size):
            checks += 1
            product_value = prod(matrix[i][j] for matrix in grams)
            expected = 1 if i == j else 0
            if product_value != expected:
                errors.append(
                    f"Hadamard product entry ({i},{j}) is {product_value}, expected {expected}"
                )
    return GramCharacterizationReport(
        not errors, tuple(grams), tuple(ranks), checks, tuple(errors[:20])
    )


def _component_extension_candidates(
    family: MVFamily, prime: int
) -> tuple[tuple[tuple[int, ...], tuple[int, ...], int, int], ...]:
    dimension = family.dimension
    residues_u = [tuple(value % prime for value in row) for row in family.u]
    residues_v = [tuple(value % prime for value in row) for row in family.v]
    vectors = tuple(cartesian_product(range(prime), repeat=dimension))
    candidates: list[tuple[tuple[int, ...], tuple[int, ...], int, int]] = []
    for new_u in vectors:
        for new_v in vectors:
            if sum(a * b for a, b in zip(new_u, new_v, strict=True)) % prime != 1:
                continue
            outgoing_zero_mask = 0
            incoming_zero_mask = 0
            valid = True
            for label, (old_u, old_v) in enumerate(zip(residues_u, residues_v, strict=True)):
                outgoing = sum(
                    a * b for a, b in zip(new_u, old_v, strict=True)
                ) % prime
                incoming = sum(
                    a * b for a, b in zip(old_u, new_v, strict=True)
                ) % prime
                if outgoing not in (0, 1) or incoming not in (0, 1):
                    valid = False
                    break
                if outgoing == 0:
                    outgoing_zero_mask |= 1 << label
                if incoming == 0:
                    incoming_zero_mask |= 1 << label
            if valid:
                candidates.append(
                    (new_u, new_v, outgoing_zero_mask, incoming_zero_mask)
                )
    return tuple(candidates)


def extend_hpr_family_one_label(family: MVFamily) -> ExtensionResult:
    """Exhaustively test whether a fixed HPR family admits one more label.

    This is not a global nonexistence test: failure proves only that the given
    witness is maximal under extension without changing its existing vectors.
    """

    base_report = verify_gram_characterization(family)
    if not base_report.ok:
        raise ValueError(f"base family is not an HPR Gram system: {base_report.errors}")
    candidates = tuple(
        _component_extension_candidates(family, prime) for prime in family.prime_factors
    )
    full_mask = (1 << family.size) - 1
    combinations_checked = 0

    def search(
        component: int,
        outgoing_mask: int,
        incoming_mask: int,
        chosen: tuple[tuple[tuple[int, ...], tuple[int, ...], int, int], ...],
    ) -> tuple[tuple[tuple[int, ...], tuple[int, ...], int, int], ...] | None:
        nonlocal combinations_checked
        if component == len(candidates):
            combinations_checked += 1
            if outgoing_mask == full_mask and incoming_mask == full_mask:
                return chosen
            return None
        for candidate in candidates[component]:
            result = search(
                component + 1,
                outgoing_mask | candidate[2],
                incoming_mask | candidate[3],
                chosen + (candidate,),
            )
            if result is not None:
                return result
        return None

    chosen = search(0, 0, 0, ())
    if chosen is None:
        return ExtensionResult(
            False, None, tuple(map(len, candidates)), combinations_checked
        )
    new_u = tuple(
        crt([chosen[index][0][coordinate] for index in range(len(chosen))], family.prime_factors)
        for coordinate in range(family.dimension)
    )
    new_v = tuple(
        crt([chosen[index][1][coordinate] for index in range(len(chosen))], family.prime_factors)
        for coordinate in range(family.dimension)
    )
    extended = MVFamily.build(
        modulus=family.modulus,
        u=family.u + (new_u,),
        v=family.v + (new_v,),
        allowed_offdiagonal=family.allowed_offdiagonal,
        diagonal_value=1,
        prime_factors=family.prime_factors,
        convention="hpr",
        metadata={**family.metadata, "extended_by_exhaustive_component_mask_search": True},
    )
    verification = verify_gram_characterization(extended)
    if not verification.ok:
        raise RuntimeError(f"extended family failed verification: {verification.errors}")
    return ExtensionResult(
        True, extended, tuple(map(len, candidates)), combinations_checked
    )


def _variable(kind: str, component: int, label: int, coordinate: int) -> str:
    return f"{kind}_{component}_{label}_{coordinate}"


def _dot_expression(
    component: int, left: int, right: int, dimension: int, prime: int
) -> str:
    terms = [
        f"(* {_variable('u', component, left, coordinate)} "
        f"{_variable('v', component, right, coordinate)})"
        for coordinate in range(dimension)
    ]
    total = terms[0] if len(terms) == 1 else f"(+ {' '.join(terms)})"
    return f"(mod {total} {prime})"


def canonical_gram_smt2(
    *,
    primes: Sequence[int],
    dimension: int,
    size: int,
    symmetry_break: bool = True,
    basis_component: int | None = None,
    other_rank_bound: int | None = None,
    other_basis_labels: Sequence[int] | None = None,
) -> tuple[str, tuple[str, ...], int]:
    """Return a complete QF_NIA factor encoding for an HPR Gram system."""

    prime_tuple = tuple(primes)
    if len(prime_tuple) < 2 or len(set(prime_tuple)) != len(prime_tuple):
        raise ValueError("at least two distinct primes are required")
    if any(prime < 3 or prime % 2 == 0 for prime in prime_tuple):
        raise ValueError("HPR requires distinct odd primes")
    if dimension < 1 or size < 1:
        raise ValueError("dimension and size must be positive")
    if basis_component is not None and not 0 <= basis_component < len(prime_tuple):
        raise ValueError("basis_component is outside the prime-factor range")
    if basis_component is not None and size < dimension:
        raise ValueError("a normalized row basis requires size >= dimension")
    if (other_rank_bound is not None or other_basis_labels is not None) and (
        basis_component is None or len(prime_tuple) != 2
    ):
        raise ValueError("other-component cases require two primes and a basis component")
    if other_rank_bound is not None and not 1 <= other_rank_bound <= dimension:
        raise ValueError("other_rank_bound is outside the dimension")
    if other_basis_labels is not None and (
        len(other_basis_labels) != dimension
        or other_basis_labels[0] != 0
        or len(set(other_basis_labels)) != dimension
        or any(not 0 <= label < size for label in other_basis_labels)
    ):
        raise ValueError("other_basis_labels must be distinct, start at zero, and have length d")
    if other_rank_bound is not None and other_basis_labels is not None:
        raise ValueError("choose a rank-bounded or a full-basis other-component case")
    other_component = None if basis_component is None else 1 - basis_component

    names = tuple(
        _variable(kind, component, label, coordinate)
        for component in range(len(prime_tuple))
        for kind in ("u", "v")
        for label in range(size)
        for coordinate in range(dimension)
    )
    lines = ["(set-logic QF_NIA)", "(set-option :produce-models true)"]
    lines.extend(f"(declare-fun {name} () Int)" for name in names)
    constraints = 0
    for component, prime in enumerate(prime_tuple):
        for kind in ("u", "v"):
            for label in range(size):
                for coordinate in range(dimension):
                    name = _variable(kind, component, label, coordinate)
                    lines.append(f"(assert (and (<= 0 {name}) (< {name} {prime})))")
                    constraints += 1

        # Since the first diagonal product is one, u_0 is nonzero.  An
        # invertible change of basis sends it to e_1 without changing the Gram
        # matrix, so these constraints lose no solutions.
        if symmetry_break:
            lines.append(f"(assert (= {_variable('u', component, 0, 0)} 1))")
            constraints += 1
            for coordinate in range(1, dimension):
                lines.append(
                    f"(assert (= {_variable('u', component, 0, coordinate)} 0))"
                )
                constraints += 1
        if symmetry_break and component == basis_component:
            # This case asserts that the chosen component has full row rank.
            # A common label permutation can put a basis containing u_0 in
            # rows 0..d-1, and an invertible basis change sends those rows to
            # the standard basis.  Searches over all possible full-rank
            # components therefore cover every solution in which some
            # component has rank d.
            for label in range(1, dimension):
                for coordinate in range(dimension):
                    value = int(label == coordinate)
                    lines.append(
                        f"(assert (= {_variable('u', component, label, coordinate)} {value}))"
                    )
                    constraints += 1
            for label in range(size):
                for coordinate in range(dimension):
                    name = _variable("v", component, label, coordinate)
                    lines.append(f"(assert (or (= {name} 0) (= {name} 1)))")
                    constraints += 1
            for left in range(dimension, size - 1):
                right = left + 1
                clauses: list[str] = []
                equal_prefix: list[str] = []
                for coordinate in range(dimension):
                    first = _variable("u", component, left, coordinate)
                    second = _variable("u", component, right, coordinate)
                    prefix = " ".join(equal_prefix)
                    clauses.append(
                        f"(and {prefix} (< {first} {second}))"
                        if prefix
                        else f"(< {first} {second})"
                    )
                    equal_prefix.append(f"(= {first} {second})")
                clauses.append(f"(and {' '.join(equal_prefix)})")
                lines.append(f"(assert (or {' '.join(clauses)}))")
                constraints += 1
        if symmetry_break and component == other_component and other_rank_bound is not None:
            for kind in ("u", "v"):
                for label in range(size):
                    for coordinate in range(other_rank_bound, dimension):
                        lines.append(
                            f"(assert (= {_variable(kind, component, label, coordinate)} 0))"
                        )
                        constraints += 1
        if symmetry_break and component == other_component and other_basis_labels is not None:
            for coordinate, label in enumerate(other_basis_labels):
                for column in range(dimension):
                    value = int(coordinate == column)
                    lines.append(
                        f"(assert (= {_variable('u', component, label, column)} {value}))"
                    )
                    constraints += 1
            for label in range(size):
                for coordinate in range(dimension):
                    name = _variable("v", component, label, coordinate)
                    lines.append(f"(assert (or (= {name} 0) (= {name} 1)))")
                    constraints += 1

        for label in range(size):
            dot = _dot_expression(component, label, label, dimension, prime)
            lines.append(f"(assert (= {dot} 1))")
            constraints += 1

    for left in range(size):
        for right in range(size):
            if left == right:
                continue
            dots = [
                _dot_expression(component, left, right, dimension, prime)
                for component, prime in enumerate(prime_tuple)
            ]
            for dot in dots:
                lines.append(f"(assert (or (= {dot} 0) (= {dot} 1)))")
                constraints += 1
            lines.append(f"(assert (or {' '.join(f'(= {dot} 0)' for dot in dots)}))")
            constraints += 1

    lines.extend(["(check-sat)", f"(get-value ({' '.join(names)}))", "(exit)"])
    return "\n".join(lines) + "\n", names, constraints


def _bit_vector_constant(value: int, width: int) -> str:
    return f"(_ bv{value} {width})"


def _bit_vector_dot_expression(
    component: int,
    left: int,
    right: int,
    dimension: int,
    prime: int,
    width: int,
) -> str:
    terms = [
        f"(bvmul {_variable('u', component, left, coordinate)} "
        f"{_variable('v', component, right, coordinate)})"
        for coordinate in range(dimension)
    ]
    total = reduce(lambda first, second: f"(bvadd {first} {second})", terms)
    return f"(bvurem {total} {_bit_vector_constant(prime, width)})"


def canonical_gram_bv_smt2(
    *,
    primes: Sequence[int],
    dimension: int,
    size: int,
    symmetry_break: bool = True,
    basis_component: int | None = None,
    other_rank_bound: int | None = None,
    other_basis_labels: Sequence[int] | None = None,
) -> tuple[str, tuple[str, ...], int]:
    """Return an exact finite bit-vector encoding of the factor search.

    Each component receives enough bits to hold the unreduced dot product, so
    bit-vector multiplication and addition do not overflow before reduction.
    Unlike the QF_NIA version, this encoding is fully bit-blastable.
    """

    prime_tuple = tuple(primes)
    if len(prime_tuple) < 2 or len(set(prime_tuple)) != len(prime_tuple):
        raise ValueError("at least two distinct primes are required")
    if any(prime < 3 or prime % 2 == 0 for prime in prime_tuple):
        raise ValueError("HPR requires distinct odd primes")
    if dimension < 1 or size < 1:
        raise ValueError("dimension and size must be positive")
    if basis_component is not None and not 0 <= basis_component < len(prime_tuple):
        raise ValueError("basis_component is outside the prime-factor range")
    if basis_component is not None and size < dimension:
        raise ValueError("a normalized row basis requires size >= dimension")
    if (other_rank_bound is not None or other_basis_labels is not None) and (
        basis_component is None or len(prime_tuple) != 2
    ):
        raise ValueError("other-component cases require two primes and a basis component")
    if other_rank_bound is not None and not 1 <= other_rank_bound <= dimension:
        raise ValueError("other_rank_bound is outside the dimension")
    if other_basis_labels is not None and (
        len(other_basis_labels) != dimension
        or other_basis_labels[0] != 0
        or len(set(other_basis_labels)) != dimension
        or any(not 0 <= label < size for label in other_basis_labels)
    ):
        raise ValueError("other_basis_labels must be distinct, start at zero, and have length d")
    if other_rank_bound is not None and other_basis_labels is not None:
        raise ValueError("choose a rank-bounded or a full-basis other-component case")
    other_component = None if basis_component is None else 1 - basis_component
    widths = tuple((dimension * (prime - 1) ** 2).bit_length() for prime in prime_tuple)
    names = tuple(
        _variable(kind, component, label, coordinate)
        for component in range(len(prime_tuple))
        for kind in ("u", "v")
        for label in range(size)
        for coordinate in range(dimension)
    )
    lines = ["(set-logic QF_BV)", "(set-option :produce-models true)"]
    for component, width in enumerate(widths):
        for kind in ("u", "v"):
            for label in range(size):
                for coordinate in range(dimension):
                    name = _variable(kind, component, label, coordinate)
                    lines.append(f"(declare-fun {name} () (_ BitVec {width}))")
    constraints = 0
    for component, (prime, width) in enumerate(zip(prime_tuple, widths, strict=True)):
        prime_constant = _bit_vector_constant(prime, width)
        for kind in ("u", "v"):
            for label in range(size):
                for coordinate in range(dimension):
                    name = _variable(kind, component, label, coordinate)
                    lines.append(f"(assert (bvult {name} {prime_constant}))")
                    constraints += 1
        if symmetry_break:
            lines.append(
                f"(assert (= {_variable('u', component, 0, 0)} "
                f"{_bit_vector_constant(1, width)}))"
            )
            constraints += 1
            for coordinate in range(1, dimension):
                lines.append(
                    f"(assert (= {_variable('u', component, 0, coordinate)} "
                    f"{_bit_vector_constant(0, width)}))"
                )
                constraints += 1
        if symmetry_break and component == basis_component:
            for label in range(1, dimension):
                for coordinate in range(dimension):
                    value = int(label == coordinate)
                    lines.append(
                        f"(assert (= {_variable('u', component, label, coordinate)} "
                        f"{_bit_vector_constant(value, width)}))"
                    )
                    constraints += 1
            zero = _bit_vector_constant(0, width)
            one = _bit_vector_constant(1, width)
            for label in range(size):
                for coordinate in range(dimension):
                    name = _variable("v", component, label, coordinate)
                    lines.append(f"(assert (or (= {name} {zero}) (= {name} {one})))")
                    constraints += 1
            for left in range(dimension, size - 1):
                right = left + 1
                clauses: list[str] = []
                equal_prefix: list[str] = []
                for coordinate in range(dimension):
                    first = _variable("u", component, left, coordinate)
                    second = _variable("u", component, right, coordinate)
                    prefix = " ".join(equal_prefix)
                    comparison = f"(bvult {first} {second})"
                    clauses.append(
                        f"(and {prefix} {comparison})" if prefix else comparison
                    )
                    equal_prefix.append(f"(= {first} {second})")
                clauses.append(f"(and {' '.join(equal_prefix)})")
                lines.append(f"(assert (or {' '.join(clauses)}))")
                constraints += 1
        if symmetry_break and component == other_component and other_rank_bound is not None:
            zero = _bit_vector_constant(0, width)
            for kind in ("u", "v"):
                for label in range(size):
                    for coordinate in range(other_rank_bound, dimension):
                        lines.append(
                            f"(assert (= {_variable(kind, component, label, coordinate)} {zero}))"
                        )
                        constraints += 1
        if symmetry_break and component == other_component and other_basis_labels is not None:
            zero = _bit_vector_constant(0, width)
            one = _bit_vector_constant(1, width)
            for coordinate, label in enumerate(other_basis_labels):
                for column in range(dimension):
                    value = one if coordinate == column else zero
                    lines.append(
                        f"(assert (= {_variable('u', component, label, column)} {value}))"
                    )
                    constraints += 1
            for label in range(size):
                for coordinate in range(dimension):
                    name = _variable("v", component, label, coordinate)
                    lines.append(f"(assert (or (= {name} {zero}) (= {name} {one})))")
                    constraints += 1
        for label in range(size):
            dot = _bit_vector_dot_expression(
                component, label, label, dimension, prime, width
            )
            lines.append(f"(assert (= {dot} {_bit_vector_constant(1, width)}))")
            constraints += 1

    for left in range(size):
        for right in range(size):
            if left == right:
                continue
            dots = [
                _bit_vector_dot_expression(
                    component,
                    left,
                    right,
                    dimension,
                    prime,
                    widths[component],
                )
                for component, prime in enumerate(prime_tuple)
            ]
            for component, dot in enumerate(dots):
                zero = _bit_vector_constant(0, widths[component])
                one = _bit_vector_constant(1, widths[component])
                lines.append(f"(assert (or (= {dot} {zero}) (= {dot} {one})))")
                constraints += 1
            zero_tests = [
                f"(= {dot} {_bit_vector_constant(0, widths[component])})"
                for component, dot in enumerate(dots)
            ]
            lines.append(f"(assert (or {' '.join(zero_tests)}))")
            constraints += 1
    lines.extend(["(check-sat)", f"(get-value ({' '.join(names)}))", "(exit)"])
    return "\n".join(lines) + "\n", names, constraints


def _parse_values(output: str, names: Sequence[str]) -> dict[str, int]:
    values: dict[str, int] = {}
    pattern = re.compile(
        r"\((u_\d+_\d+_\d+|v_\d+_\d+_\d+)\s+"
        r"(#b[01]+|#x[0-9A-Fa-f]+|-?\d+|\(_\s+bv\d+\s+\d+\))\)"
    )
    for name, raw_value in pattern.findall(output):
        if raw_value.startswith("#b"):
            values[name] = int(raw_value[2:], 2)
        elif raw_value.startswith("#x"):
            values[name] = int(raw_value[2:], 16)
        elif raw_value.startswith("(_"):
            values[name] = int(re.search(r"bv(\d+)", raw_value).group(1))
        else:
            values[name] = int(raw_value)
    missing = set(names) - values.keys()
    if missing:
        raise RuntimeError(f"solver model omitted {len(missing)} variables")
    return values


def _family_from_values(
    *, values: dict[str, int], primes: tuple[int, ...], dimension: int, size: int
) -> MVFamily:
    component_u = [
        [
            tuple(
                values[_variable("u", component, label, coordinate)]
                for coordinate in range(dimension)
            )
            for label in range(size)
        ]
        for component in range(len(primes))
    ]
    component_v = [
        [
            tuple(
                values[_variable("v", component, label, coordinate)]
                for coordinate in range(dimension)
            )
            for label in range(size)
        ]
        for component in range(len(primes))
    ]
    u = [
        tuple(
            crt(
                [component_u[index][label][coordinate] for index in range(len(primes))],
                primes,
            )
            for coordinate in range(dimension)
        )
        for label in range(size)
    ]
    v = [
        tuple(
            crt(
                [component_v[index][label][coordinate] for index in range(len(primes))],
                primes,
            )
            for coordinate in range(dimension)
        )
        for label in range(size)
    ]
    modulus = prod(primes)
    return MVFamily.build(
        modulus=modulus,
        u=u,
        v=v,
        allowed_offdiagonal=canonical_hpr_offdiagonal_set(primes),
        diagonal_value=1,
        prime_factors=primes,
        convention="hpr",
        metadata={
            "construction": "direct SMT factorization of binary CRT Gram matrices",
            "size": size,
            "dimension": dimension,
        },
    )


def run_canonical_gram_smt(
    *,
    primes: Sequence[int],
    dimension: int,
    size: int,
    timeout_seconds: int = 120,
    symmetry_break: bool = True,
    encoding: str = "bv",
    solver_seed: int = 0,
    basis_component: int | None = None,
    other_rank_bound: int | None = None,
    other_basis_labels: Sequence[int] | None = None,
) -> GramSMTResult:
    """Search exactly for a bounded HPR family and verify any witness."""

    solver_path = shutil.which("z3")
    if solver_path is None:
        raise RuntimeError("z3 executable not found")
    prime_tuple = tuple(primes)
    encoders: dict[str, Callable[..., tuple[str, tuple[str, ...], int]]] = {
        "bv": canonical_gram_bv_smt2,
        "nia": canonical_gram_smt2,
    }
    if encoding not in encoders:
        raise ValueError("encoding must be 'bv' or 'nia'")
    smt_encoding, names, constraints = encoders[encoding](
        primes=prime_tuple,
        dimension=dimension,
        size=size,
        symmetry_break=symmetry_break,
        basis_component=basis_component,
        other_rank_bound=other_rank_bound,
        other_basis_labels=other_basis_labels,
    )
    digest = sha256(smt_encoding.encode("utf-8")).hexdigest()
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix="treeeval-gram-smt-") as directory:
        path = Path(directory) / "instance.smt2"
        path.write_text(smt_encoding, encoding="utf-8")
        process = subprocess.run(
            [
                solver_path,
                f"-T:{timeout_seconds}",
                f"sat.random_seed={solver_seed}",
                f"smt.random_seed={solver_seed}",
                str(path),
            ],
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout_seconds + 10,
        )
    wall_seconds = time.monotonic() - started
    output = process.stdout.strip()
    first_line = output.splitlines()[0] if output else f"error(exit={process.returncode})"
    family = None
    verification = None
    if first_line == "sat":
        values = _parse_values(output, names)
        family = _family_from_values(
            values=values, primes=prime_tuple, dimension=dimension, size=size
        )
        verification = verify_gram_characterization(family)
        if not verification.ok:
            raise RuntimeError(
                f"solver witness failed Gram verification: {verification.errors}"
            )
    version = subprocess.run(
        [solver_path, "-version"], check=False, capture_output=True, text=True
    ).stdout.strip()
    return GramSMTResult(
        status=first_line,
        size=size,
        dimension=dimension,
        primes=prime_tuple,
        encoding_sha256=digest,
        solver=version,
        wall_seconds=wall_seconds,
        variables=len(names),
        constraints=constraints,
        solver_seed=solver_seed,
        family=family,
        verification=verification,
    )
