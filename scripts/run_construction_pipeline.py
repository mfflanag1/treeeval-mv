#!/usr/bin/env python3
"""Reproduce the construction -> changed type -> verification pipeline."""

from __future__ import annotations

import json

from treeeval_mv.constructions import coordinate_equality_hpr_family
from treeeval_mv.factored_vector import (
    FactoredEqualityVectors,
    Overlay,
    apply_overlay,
    query_overlay_direct,
    query_overlay_primewise,
    recursive_type_accounting,
)
from treeeval_mv.gram_search import verify_gram_characterization
from treeeval_mv.independent_verifier import verify_mv_family_independently
from treeeval_mv.structured_kernels import (
    kernel_pair_report,
    kernel_pair_to_hpr_family,
    symplectic_oriented_kernels,
)
from treeeval_mv.verifier import verify_mv_family


def main() -> None:
    equality = coordinate_equality_hpr_family(
        alphabet_size=3, primes=(3, 5, 7)
    )
    equality_checks = {
        "primary_mv": verify_mv_family(equality).ok,
        "independent_mv": verify_mv_family_independently(equality).ok,
        "gram": verify_gram_characterization(equality).ok,
    }

    first, second = symplectic_oriented_kernels(1)
    structured_report = kernel_pair_report(first, second, (3, 5))
    structured = kernel_pair_to_hpr_family(first, second, (3, 5))
    structured_checks = {
        "primary_mv": verify_mv_family(structured).ok,
        "independent_mv": verify_mv_family_independently(structured).ok,
        "gram": verify_gram_characterization(structured).ok,
    }

    vectors = FactoredEqualityVectors(3, (3, 5, 7))
    base = (17, 83, 41)
    overlay_checks = 0
    reversal_checks = 0
    for left in range(vectors.capacity):
        overlay = Overlay(left, 19)
        for right in range(vectors.capacity):
            if query_overlay_direct(vectors, base, overlay, right) != query_overlay_primewise(
                vectors, base, overlay, right
            ):
                raise AssertionError("independent overlay queries disagree")
            overlay_checks += 1
        changed = apply_overlay(vectors, base, overlay)
        if apply_overlay(vectors, changed, Overlay(left, -19)) != base:
            raise AssertionError("balanced overlay failed to reverse")
        reversal_checks += 1

    asymptotic_samples = []
    for scale in (4, 8, 16, 32):
        accounting = recursive_type_accounting(
            height=scale,
            label_bits=scale,
            prime_count=scale,
            calls_per_level=4**scale,
        )
        asymptotic_samples.append(
            {
                "ell_equals_h": scale,
                "branch_exponent": accounting.branch_exponent,
                "input_log_scale": 2 * scale,
                "polynomial_degree_required": accounting.branch_exponent // (2 * scale),
            }
        )

    output = {
        "stage_1": {
            "coordinate_equality": {
                "N": equality.size,
                "d": equality.dimension,
                "m": equality.modulus,
                "checks": equality_checks,
            },
            "oriented_symplectic": {
                "N": structured_report.size,
                "d": structured_report.dimension,
                "ranks": structured_report.ranks,
                "checks": structured_checks,
                "beats_equality_baseline": False,
            },
        },
        "stage_2": {
            "implicit_coordinate_matches": all(
                vectors.materialize(label) == equality.u[label]
                for label in range(vectors.capacity)
            ),
            "independent_overlay_checks": overlay_checks,
            "exact_reversal_checks": reversal_checks,
        },
        "stage_3": {
            "exact_algebra": True,
            "two_independent_verifiers": True,
            "exact_reversal": True,
            "polynomial_time_for_growing_t": False,
            "reason": "required polynomial degree grows with ell",
            "samples": asymptotic_samples,
        },
    }
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
