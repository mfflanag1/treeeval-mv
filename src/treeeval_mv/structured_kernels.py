"""Construction-first search over characteristic-dependent binary kernels.

Every binary matrix pair ``G,H`` with diagonal one, ``G circ H = I``, and
small ranks in two different characteristics yields a canonical HPR matching
vector family.  This module keeps the search at the Gram level, then performs
an exact rank factorization and CRT lift before a candidate is accepted.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from math import prod
from typing import Sequence

from .gram_search import Matrix, matrix_rank_mod
from .model import MVFamily, canonical_hpr_offdiagonal_set, crt


@dataclass(frozen=True)
class KernelPairReport:
    size: int
    primes: tuple[int, int]
    ranks: tuple[int, int]
    dimension: int
    hadamard_identity: bool
    binary_diagonal_one: bool


def _inverse_mod(matrix: Sequence[Sequence[int]], prime: int) -> list[list[int]]:
    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise ValueError("matrix must be square")
    work = [
        [value % prime for value in row]
        + [int(i == j) for j in range(size)]
        for i, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = next((r for r in range(column, size) if work[r][column]), None)
        if pivot is None:
            raise ValueError("matrix is singular")
        work[column], work[pivot] = work[pivot], work[column]
        inverse = pow(work[column][column], -1, prime)
        work[column] = [(value * inverse) % prime for value in work[column]]
        for row in range(size):
            if row == column or not work[row][column]:
                continue
            multiplier = work[row][column]
            work[row] = [
                (left - multiplier * right) % prime
                for left, right in zip(work[row], work[column], strict=True)
            ]
    return [row[size:] for row in work]


def factor_matrix_mod(matrix: Matrix, prime: int) -> tuple[Matrix, Matrix]:
    """Return ``U,V`` with ``matrix[i][j]=<U[i],V[j]> mod prime``."""

    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise ValueError("only square Gram matrices are supported")
    rank = matrix_rank_mod(matrix, prime)
    if rank == 0:
        zero = tuple(() for _ in range(size))
        return zero, zero

    pivot_columns: list[int] = []
    for column in range(size):
        candidate = pivot_columns + [column]
        columns_as_rows = tuple(
            tuple(matrix[row][index] % prime for row in range(size))
            for index in candidate
        )
        if matrix_rank_mod(columns_as_rows, prime) > len(pivot_columns):
            pivot_columns.append(column)
        if len(pivot_columns) == rank:
            break
    left = tuple(
        tuple(matrix[row][column] % prime for column in pivot_columns)
        for row in range(size)
    )

    pivot_rows: list[int] = []
    for row in range(size):
        candidate = pivot_rows + [row]
        if matrix_rank_mod(tuple(left[index] for index in candidate), prime) > len(
            pivot_rows
        ):
            pivot_rows.append(row)
        if len(pivot_rows) == rank:
            break
    square = tuple(left[row] for row in pivot_rows)
    inverse = _inverse_mod(square, prime)
    right_rows: list[tuple[int, ...]] = []
    for column in range(size):
        restricted = [matrix[row][column] % prime for row in pivot_rows]
        coordinates = tuple(
            sum(inverse[i][j] * restricted[j] for j in range(rank)) % prime
            for i in range(rank)
        )
        right_rows.append(coordinates)

    right = tuple(right_rows)
    for i in range(size):
        for j in range(size):
            value = sum(a * b for a, b in zip(left[i], right[j], strict=True)) % prime
            if value != matrix[i][j] % prime:
                raise AssertionError("rank factorization failed its exact reconstruction check")
    return left, right


def kernel_pair_report(
    first: Matrix, second: Matrix, primes: tuple[int, int]
) -> KernelPairReport:
    if len(first) != len(second) or any(len(row) != len(first) for row in first + second):
        raise ValueError("kernel matrices must be square and have equal size")
    size = len(first)
    binary_diagonal_one = all(
        first[i][j] in (0, 1)
        and second[i][j] in (0, 1)
        and (i != j or (first[i][j] == second[i][j] == 1))
        for i in range(size)
        for j in range(size)
    )
    hadamard_identity = all(
        first[i][j] * second[i][j] == int(i == j)
        for i in range(size)
        for j in range(size)
    )
    ranks = (
        matrix_rank_mod(first, primes[0]),
        matrix_rank_mod(second, primes[1]),
    )
    return KernelPairReport(
        size, primes, ranks, max(ranks), hadamard_identity, binary_diagonal_one
    )


def kernel_pair_to_hpr_family(
    first: Matrix, second: Matrix, primes: tuple[int, int]
) -> MVFamily:
    report = kernel_pair_report(first, second, primes)
    if not report.binary_diagonal_one or not report.hadamard_identity:
        raise ValueError("kernels do not form a canonical binary Gram pair")
    factors = (
        factor_matrix_mod(first, primes[0]),
        factor_matrix_mod(second, primes[1]),
    )
    dimension = report.dimension

    def lift(side: int, label: int) -> tuple[int, ...]:
        return tuple(
            crt(
                [
                    factors[component][side][label][coordinate]
                    if coordinate < len(factors[component][side][label])
                    else 0
                    for component in range(2)
                ],
                primes,
            )
            for coordinate in range(dimension)
        )

    modulus = prod(primes)
    return MVFamily.build(
        modulus=modulus,
        u=(lift(0, label) for label in range(report.size)),
        v=(lift(1, label) for label in range(report.size)),
        allowed_offdiagonal=canonical_hpr_offdiagonal_set(primes),
        diagonal_value=1,
        prime_factors=primes,
        convention="hpr",
        metadata={
            "construction": "exact factorization of structured binary kernels",
            "component_ranks": report.ranks,
        },
    )


def symplectic_oriented_kernels(pair_count: int) -> tuple[Matrix, Matrix]:
    """The zero-form kernel and one oriented nonzero-form kernel over F_3."""

    if pair_count < 1:
        raise ValueError("pair_count must be positive")
    vectors = tuple(product(range(3), repeat=2 * pair_count))

    def form(left: tuple[int, ...], right: tuple[int, ...]) -> int:
        return sum(
            left[2 * i] * right[2 * i + 1]
            - left[2 * i + 1] * right[2 * i]
            for i in range(pair_count)
        ) % 3

    first = tuple(
        tuple(int(form(left, right) == 0) for right in vectors) for left in vectors
    )
    second = tuple(
        tuple(int(i == j or form(left, right) == 1) for j, right in enumerate(vectors))
        for i, left in enumerate(vectors)
    )
    return first, second
