#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

from treeeval_mv.bit_tape_packing import exhaust_bit_tape_transitions
from treeeval_mv.constructions import (
    dgy_grolmusz_family,
    dgy_grolmusz_polynomial_family,
)
from treeeval_mv.hpr_selector import verify_hpr_selector_identity
from treeeval_mv.gks_cir_calibration import verify_mod6_derivative_cir
from treeeval_mv.gram_search import equality_product_hpr_family
from treeeval_mv.independent_verifier import verify_mv_family_independently
from treeeval_mv.model import flip_customary_to_hpr
from treeeval_mv.packed_catalyst import exhaust_packing_state_map, packed_bit_length
from treeeval_mv.projected_selector import (
    has_cross_characteristic_fixed_weight,
    verify_projected_selector_identity,
)
from treeeval_mv.succinct_profile_barrier import (
    dynamic_profile_barrier,
    zero_base_scan_report,
)
from treeeval_mv.selector_hitting_set import minimum_selector_hitting_set
from treeeval_mv.verifier import verify_mv_family, verify_polynomial_matching_family


ROOT = Path(__file__).resolve().parents[1]


def report_dict(report: object) -> dict[str, object]:
    return {
        "ok": getattr(report, "ok"),
        "checks": getattr(report, "checks"),
        "errors": list(getattr(report, "errors")),
    }


def family_dict(family: object) -> dict[str, object]:
    return {
        "modulus": getattr(family, "modulus"),
        "size": getattr(family, "size"),
        "dimension": getattr(family, "dimension"),
        "prime_factors": list(getattr(family, "prime_factors")),
        "allowed_offdiagonal": sorted(getattr(family, "allowed_offdiagonal")),
        "u": [list(row) for row in getattr(family, "u")],
        "v": [list(row) for row in getattr(family, "v")],
        "metadata": dict(getattr(family, "metadata")),
    }


def packing_report_dict(report: object) -> dict[str, object]:
    payload = {
        "ok": getattr(report, "ok"),
        "initial_states": getattr(report, "initial_states"),
        "errors": list(getattr(report, "errors")),
    }
    if hasattr(report, "update_checks"):
        payload["update_checks"] = getattr(report, "update_checks")
    if hasattr(report, "transitions"):
        payload["transitions"] = getattr(report, "transitions")
    return payload


def main() -> None:
    polynomial6 = dgy_grolmusz_polynomial_family(
        primes=(2, 3), exponents=(1, 1), weight=2, universe_size=4
    )
    family6 = dgy_grolmusz_family(
        primes=(2, 3), exponents=(1, 1), weight=2, universe_size=4
    )
    family15 = dgy_grolmusz_family(
        primes=(3, 5), exponents=(1, 1), weight=2, universe_size=4
    )
    hpr15 = flip_customary_to_hpr(family15)
    derivative_cir = verify_mod6_derivative_cir(family6)
    projected_selector = verify_projected_selector_identity(hpr15)
    profile_family = equality_product_hpr_family(primes=(3, 5), alphabet_size=2)
    profile_barrier = dynamic_profile_barrier(profile_family)
    hpr15_profile_barrier = dynamic_profile_barrier(hpr15)
    zero_scan = zero_base_scan_report(profile_family)
    selector_hitting_sets = tuple(
        minimum_selector_hitting_set(prime) for prime in (3, 5, 7, 11, 13)
    )
    output = {
        "mod_6_dgy_grolmusz_calibration": {
            "family": family_dict(family6),
            "polynomial_metadata": dict(polynomial6.metadata),
            "polynomial_matching_verifier": report_dict(
                verify_polynomial_matching_family(polynomial6)
            ),
            "mv_verifier_primary": report_dict(verify_mv_family(family6)),
            "mv_verifier_independent_crt": report_dict(
                verify_mv_family_independently(family6)
            ),
            "hpr_selector": {
                "status": "not-applicable",
                "reason": "HPR assumes odd primes; selector coefficients ±2 must be units",
            },
        },
        "mod_15_hpr_selector_companion": {
            "family": family_dict(hpr15),
            "mv_verifier_primary": report_dict(verify_mv_family(hpr15)),
            "mv_verifier_independent_crt": report_dict(
                verify_mv_family_independently(hpr15)
            ),
            "hpr_lemmas_3_7_to_3_10_selector": report_dict(
                verify_hpr_selector_identity(hpr15)
            ),
            "characteristic_local_projected_selector": {
                "ok": projected_selector.ok,
                "brute_force_checks": projected_selector.brute_force_checks,
                "factor_checks": projected_selector.factor_checks,
                "errors": list(projected_selector.errors),
            },
            "cross_characteristic_fixed_weight": {
                "F3_to_F5_exists": has_cross_characteristic_fixed_weight(3, 5),
                "F5_to_F3_exists": has_cross_characteristic_fixed_weight(5, 3),
                "expected": False,
            },
        },
        "joint_packed_catalyst_calibration": {
            "modulus": 3,
            "coordinate_count": 4,
            "catalytic_bits": packed_bit_length(3, 4),
            "quotient_free_bits": 1,
            "exhaustive_report": packing_report_dict(
                exhaust_packing_state_map(modulus=3, coordinate_count=4)
            ),
            "bit_level_transition_report": packing_report_dict(
                exhaust_bit_tape_transitions(modulus=3, coordinate_count=4)
            ),
        },
        "mod_6_derivative_pair_kernel_calibration": {
            "model": "GKS multiplicity-two PIR pair kernel, GF(4), MV modulus 6",
            "ok": derivative_cir.ok,
            "single_kernel_checks": derivative_cir.single_kernel_checks,
            "adversarial_direction_mask_checks": derivative_cir.derivative_mask_checks,
            "pair_kernel_checks": derivative_cir.pair_kernel_checks,
            "scope": "algebra only: every source/target linear kernel coefficient under two exponent masks and two direction masks; excludes recursive field-to-exponent conversion",
            "errors": list(derivative_cir.errors),
        },
        "dynamic_mv_profile_barrier_calibration": {
            "family": "coordinate-equality product over Z_15, alphabet size 2",
            "component_ranks": list(profile_barrier.component_ranks),
            "component_profile_counts": list(
                profile_barrier.component_profile_counts
            ),
            "total_profile_count": profile_barrier.total_profile_count,
            "minimum_bits": profile_barrier.minimum_bits,
            "exhaustive_components": list(profile_barrier.exhaustive_components),
            "zero_base_scan": {
                "active_candidate_pairs": zero_scan.active_candidate_pairs,
                "update_scalar": zero_scan.update_scalar,
                "scalar_is_unit": zero_scan.scalar_is_unit,
            },
            "dgy_hpr_Z15_companion": {
                "component_ranks": list(hpr15_profile_barrier.component_ranks),
                "component_profile_counts": list(
                    hpr15_profile_barrier.component_profile_counts
                ),
                "total_profile_count": hpr15_profile_barrier.total_profile_count,
                "minimum_bits": hpr15_profile_barrier.minimum_bits,
                "exhaustive_components": list(
                    hpr15_profile_barrier.exhaustive_components
                ),
            },
        },
        "selector_exponent_hitting_set_calibration": {
            str(report.prime): {
                "minimum_size": report.minimum_size,
                "witness": list(report.witness),
                "base_pairs": report.base_pairs,
                "maximum_single_exponent_coverage": (
                    report.maximum_single_exponent_coverage
                ),
                "counting_lower_bound": report.counting_lower_bound,
            }
            for report in selector_hitting_sets
        },
    }
    target = ROOT / "results" / "calibration.json"
    target.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(target)


if __name__ == "__main__":
    main()
