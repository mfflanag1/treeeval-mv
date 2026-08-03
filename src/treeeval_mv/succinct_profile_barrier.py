"""Information lower bound for exact dynamic matching-vector sketches.

The bound applies to a context-independent state representation that starts at
zero, supports arbitrary matching-vector updates, and answers every Gram
query.  It does not rule out trace-specific recomputation or a representation
that exploits the control state of the TreeEval recursion.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from math import prod

from .gram_search import matrix_rank_mod, verify_gram_characterization
from .model import MVFamily


@dataclass(frozen=True)
class DynamicProfileBarrierReport:
    primes: tuple[int, ...]
    component_ranks: tuple[int, ...]
    component_profile_counts: tuple[int, ...]
    total_profile_count: int
    minimum_bits: int
    exhaustive_components: tuple[bool, ...]


@dataclass(frozen=True)
class ZeroBaseScanReport:
    family_size: int
    active_candidate_pairs: int
    update_scalar: int
    scalar_is_unit: bool


def _profile_count_by_enumeration(
    gram: tuple[tuple[int, ...], ...], prime: int
) -> int:
    """Count profiles ``c^T G`` by exhaustive coefficient enumeration."""

    profiles = {
        tuple(
            sum(coefficients[row] * gram[row][column] for row in range(len(gram)))
            % prime
            for column in range(len(gram))
        )
        for coefficients in product(range(prime), repeat=len(gram))
    }
    return len(profiles)


def dynamic_profile_barrier(
    family: MVFamily, *, exhaustive_coefficient_limit: int = 1_000_000
) -> DynamicProfileBarrierReport:
    """Return the exact state-count lower bound for a dynamic MV sketch.

    Starting from zero, arbitrary updates by the rows of ``U`` reach every
    vector in their span.  The vector of all answers against rows of ``V`` is
    ``c^T G``.  Over component field ``F_p`` these answer profiles form the
    row space of ``G`` and therefore number ``p**rank(G)``.  CRT lets the
    update coefficients be chosen independently in every component, so the
    joint number of distinguishable states is the product of those counts.

    No linearity of the proposed encoding is assumed: two different answer
    profiles must occupy different states if all queries are answered exactly.
    """

    characterization = verify_gram_characterization(family)
    if not characterization.ok:
        raise ValueError(
            "the dynamic-profile barrier requires an HPR canonical family: "
            + "; ".join(characterization.errors)
        )

    ranks: list[int] = []
    counts: list[int] = []
    exhaustive: list[bool] = []
    for prime, gram in zip(
        family.prime_factors, characterization.component_grams, strict=True
    ):
        rank = matrix_rank_mod(gram, prime)
        expected_count = prime**rank
        coefficient_count = prime**family.size
        if coefficient_count <= exhaustive_coefficient_limit:
            enumerated_count = _profile_count_by_enumeration(gram, prime)
            if enumerated_count != expected_count:
                raise AssertionError(
                    f"mod {prime}: enumerated {enumerated_count} profiles, "
                    f"rank predicts {expected_count}"
                )
            exhaustive.append(True)
        else:
            exhaustive.append(False)
        ranks.append(rank)
        counts.append(expected_count)

    total = prod(counts)
    return DynamicProfileBarrierReport(
        primes=family.prime_factors,
        component_ranks=tuple(ranks),
        component_profile_counts=tuple(counts),
        total_profile_count=total,
        minimum_bits=(total - 1).bit_length(),
        exhaustive_components=tuple(exhaustive),
    )


def zero_base_scan_report(family: MVFamily) -> ZeroBaseScanReport:
    """Describe HPR's all-zero-shift scan when both logical masks are zero.

    At every component the selector polynomial is ``X-1``, so the selected
    monomial has ``(alpha,beta)=(-1,0)``.  Every candidate pair has current
    product zero and is therefore active.  Because the TreeEval truth table is
    arbitrary, this prefix exposes an arbitrary sequence of output-vector
    labels to the evolving register.
    """

    characterization = verify_gram_characterization(family)
    if not characterization.ok:
        raise ValueError(
            "the zero-base scan report requires an HPR canonical family: "
            + "; ".join(characterization.errors)
        )
    scalar = pow((-1) ** len(family.prime_factors), -1, family.modulus)
    return ZeroBaseScanReport(
        family_size=family.size,
        active_candidate_pairs=family.size * family.size,
        update_scalar=scalar,
        scalar_is_unit=True,
    )
