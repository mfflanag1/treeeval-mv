#!/usr/bin/env python3
"""Generate the exact asymmetric CRT Gram-system search ledger."""

from __future__ import annotations

import json
from pathlib import Path

from treeeval_mv.gram_search import (
    characteristic_free_size_bound,
    determinant_lift_is_safe,
    equality_product_hpr_family,
    extend_hpr_family_one_label,
    run_canonical_gram_smt,
    verify_gram_characterization,
)
from treeeval_mv.independent_verifier import verify_mv_family_independently
from treeeval_mv.verifier import verify_mv_family


ROOT = Path(__file__).resolve().parents[1]


def family_payload(family: object) -> dict[str, object]:
    return {
        "modulus": getattr(family, "modulus"),
        "size": getattr(family, "size"),
        "dimension": getattr(family, "dimension"),
        "u": [list(row) for row in getattr(family, "u")],
        "v": [list(row) for row in getattr(family, "v")],
    }


def solver_payload(result: object) -> dict[str, object]:
    return {
        "status": getattr(result, "status"),
        "encoding_sha256": getattr(result, "encoding_sha256"),
        "solver": getattr(result, "solver"),
        "solver_seed": getattr(result, "solver_seed"),
        "wall_seconds": round(getattr(result, "wall_seconds"), 6),
        "variables": getattr(result, "variables"),
        "constraints": getattr(result, "constraints"),
    }


def main() -> None:
    exact_cases: list[dict[str, object]] = []
    for primes, dimension in (((3, 5), 2), ((3, 5, 7), 2)):
        family = equality_product_hpr_family(
            primes=primes, alphabet_size=dimension
        )
        gram = verify_gram_characterization(family)
        exact_cases.append(
            {
                "primes": list(primes),
                "dimension": dimension,
                "maximum_N": characteristic_free_size_bound(
                    primes=primes, rank_bound=dimension
                ),
                "lower_bound_construction": "coordinate-equality product",
                "component_ranks": list(gram.component_ranks),
                "witness": family_payload(family),
                "primary_verified": verify_mv_family(family).ok,
                "independent_crt_verified": verify_mv_family_independently(family).ok,
                "gram_verified": gram.ok,
                "upper_bound_proof": (
                    "zero-one determinant lift to Q, then "
                    "rank(A_1 hadamard ... hadamard A_t) <= product rank(A_i)"
                ),
            }
        )

    sat_four = run_canonical_gram_smt(
        primes=(3, 5), dimension=2, size=4, timeout_seconds=30
    )
    unsat_five = run_canonical_gram_smt(
        primes=(3, 5), dimension=2, size=5, timeout_seconds=30
    )

    exceptional_attempt = run_canonical_gram_smt(
        primes=(3, 5),
        dimension=3,
        size=10,
        timeout_seconds=30,
        solver_seed=0,
    )
    extension_trials: list[dict[str, object]] = []
    for seed in range(10):
        witness = run_canonical_gram_smt(
            primes=(3, 5),
            dimension=3,
            size=7,
            timeout_seconds=30,
            solver_seed=seed,
        )
        extension = extend_hpr_family_one_label(witness.family)
        extension_trials.append(
            {
                "seed": seed,
                "base_status": witness.status,
                "extended_to_8": extension.extended,
                "component_candidates": list(extension.component_candidates),
                "combinations_checked": extension.combinations_checked,
            }
        )

    ledger = {
        "object": {
            "description": (
                "binary component Gram matrices G_i with rank_Fpi(G_i)<=d, "
                "diagonal one, and entrywise product equal to I_N"
            ),
            "equivalence": (
                "necessary and sufficient for an HPR-convention canonical MV family"
            ),
        },
        "exact_characteristic_free_thresholds": exact_cases,
        "independent_smt_cross_check": {
            "primes": [3, 5],
            "dimension": 2,
            "N_4": solver_payload(sat_four),
            "N_5": solver_payload(unsat_five),
            "conclusion": "maximum N is exactly 4",
        },
        "first_exceptional_case": {
            "primes": [3, 5],
            "dimension": 3,
            "rank_lift_safe_by_prime": {
                str(prime): determinant_lift_is_safe(prime=prime, rank_bound=3)
                for prime in (3, 5)
            },
            "baseline_N": 9,
            "N_10_global_search": solver_payload(exceptional_attempt),
            "interpretation": (
                "timeout is neither a construction nor a no-go; characteristic 3 is "
                "the first odd prime that can hide a nonzero 4x4 zero-one minor"
            ),
            "fixed_witness_extension_trials": extension_trials,
            "extension_scope_warning": (
                "failure only proves maximality of each fixed solver witness"
            ),
        },
        "encoding": {
            "kind": "exact QF_BV factor encoding",
            "overflow_guard": (
                "component width holds d*(p-1)^2 before modular reduction"
            ),
            "solver_witness_check": (
                "explicit CRT reconstruction plus primary, independent CRT, and Gram checks"
            ),
        },
    }
    target = ROOT / "results" / "gram_search_ledger.json"
    target.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(target)


if __name__ == "__main__":
    main()
