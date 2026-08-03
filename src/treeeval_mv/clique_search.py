"""Exact bounded MV search through maximum clique.

For fixed (m,d,S), each vertex is a diagonal-valid pair (u,v).  Two vertices
are adjacent exactly when both cross inner products lie in S.  Thus MV families
are cliques, with no relaxation or arithmetic encoding gap.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, product
from typing import Iterable, Sequence

from .model import MVFamily, Vector


@dataclass(frozen=True)
class CliqueResult:
    maximum_size: int
    clique: tuple[int, ...]
    vertices: int
    edges: int
    search_nodes: int


def _dot(left: Vector, right: Vector, modulus: int) -> int:
    return sum(a * b for a, b in zip(left, right, strict=True)) % modulus


def build_mv_compatibility_graph(
    *, modulus: int, dimension: int, allowed_offdiagonal: Iterable[int], diagonal_value: int = 0
) -> tuple[tuple[tuple[Vector, Vector], ...], tuple[int, ...]]:
    allowed = frozenset(value % modulus for value in allowed_offdiagonal)
    vectors = tuple(product(range(modulus), repeat=dimension))
    dot_table = tuple(
        tuple(_dot(left, right, modulus) for right in vectors) for left in vectors
    )
    vertices = tuple(
        (left, right)
        for left_index, left in enumerate(vectors)
        for right_index, right in enumerate(vectors)
        if dot_table[left_index][right_index] == diagonal_value % modulus
    )
    vector_index = {vector: index for index, vector in enumerate(vectors)}
    adjacency: list[int] = []
    for vertex_index, (left, right) in enumerate(vertices):
        left_index = vector_index[left]
        right_index = vector_index[right]
        neighbors = 0
        for other_index, (other_left, other_right) in enumerate(vertices):
            if vertex_index == other_index:
                continue
            if (
                dot_table[left_index][vector_index[other_right]] in allowed
                and dot_table[vector_index[other_left]][right_index] in allowed
            ):
                neighbors |= 1 << other_index
        adjacency.append(neighbors)
    return vertices, tuple(adjacency)


def maximum_clique(adjacency: Sequence[int]) -> CliqueResult:
    """Tomita-style branch-and-bound with an exact greedy-color upper bound."""

    best: list[int] = []
    search_nodes = 0

    def color_sort(candidates: int) -> tuple[list[int], list[int]]:
        order: list[int] = []
        bounds: list[int] = []
        color = 0
        remaining = candidates
        while remaining:
            color += 1
            independent = remaining
            while independent:
                bit = independent & -independent
                vertex = bit.bit_length() - 1
                order.append(vertex)
                bounds.append(color)
                remaining ^= bit
                independent ^= bit
                independent &= ~adjacency[vertex]
        return order, bounds

    def expand(candidates: int, current: list[int]) -> None:
        nonlocal best, search_nodes
        search_nodes += 1
        order, color_bounds = color_sort(candidates)
        for index in range(len(order) - 1, -1, -1):
            if len(current) + color_bounds[index] <= len(best):
                return
            vertex = order[index]
            bit = 1 << vertex
            if not candidates & bit:
                continue
            new_candidates = candidates & adjacency[vertex]
            if new_candidates:
                expand(new_candidates, current + [vertex])
            elif len(current) + 1 > len(best):
                best = current + [vertex]
            candidates ^= bit

    expand((1 << len(adjacency)) - 1, [])
    edges = sum(neighbors.bit_count() for neighbors in adjacency) // 2
    return CliqueResult(len(best), tuple(best), len(adjacency), edges, search_nodes)


def exact_mv_search(
    *,
    modulus: int,
    dimension: int,
    allowed_offdiagonal: Iterable[int],
    diagonal_value: int = 0,
    prime_factors: Sequence[int] = (),
) -> tuple[MVFamily, CliqueResult, tuple[int, ...]]:
    allowed = frozenset(allowed_offdiagonal)
    vertices, adjacency = build_mv_compatibility_graph(
        modulus=modulus,
        dimension=dimension,
        allowed_offdiagonal=allowed,
        diagonal_value=diagonal_value,
    )
    result = maximum_clique(adjacency)
    pairs = [vertices[index] for index in result.clique]
    family = MVFamily.build(
        modulus=modulus,
        u=[left for left, _ in pairs],
        v=[right for _, right in pairs],
        allowed_offdiagonal=allowed,
        diagonal_value=diagonal_value,
        prime_factors=prime_factors,
        metadata={
            "construction": "exact maximum clique in the full diagonal-valid pair graph",
            "bounded_search": {"modulus": modulus, "dimension": dimension},
        },
    )
    return family, result, adjacency


def structured_set_system_search(
    *, modulus: int, universe_size: int, weight: int, allowed_intersections: Iterable[int]
) -> tuple[tuple[frozenset[int], ...], CliqueResult]:
    """Exhaust the constant-weight restricted-intersection construction class."""

    allowed = frozenset(value % modulus for value in allowed_intersections)
    if weight % modulus != 0:
        raise ValueError("set size must be zero modulo the modulus for incidence-vector diagonal zero")
    sets = tuple(frozenset(values) for values in combinations(range(universe_size), weight))
    adjacency: list[int] = []
    for index, left in enumerate(sets):
        neighbors = 0
        for other_index, right in enumerate(sets):
            if index != other_index and len(left & right) % modulus in allowed:
                neighbors |= 1 << other_index
        adjacency.append(neighbors)
    result = maximum_clique(adjacency)
    return tuple(sets[index] for index in result.clique), result
