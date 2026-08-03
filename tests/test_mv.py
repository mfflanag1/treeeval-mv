from __future__ import annotations

import unittest

from treeeval_mv.bit_tape_packing import exhaust_bit_tape_transitions
from treeeval_mv.branch_aware_barrier import (
    branch_aware_exact_barrier,
    multiplicative_unit_subgroups,
    unit_edge_potential,
    verify_unit_edge_potential,
)
from treeeval_mv.clique_search import exact_mv_search, structured_set_system_search
from treeeval_mv.constructions import (
    coordinate_equality_hpr_family,
    dgy_grolmusz_family,
    dgy_grolmusz_polynomial_family,
    weight_congruence_test,
)
from treeeval_mv.hpr_selector import verify_hpr_selector_identity
from treeeval_mv.holonomy_barrier import (
    minimum_reversible_holonomy_states,
    verify_holonomy_bound,
)
from treeeval_mv.gks_cir_calibration import (
    F4_ONE,
    f4_mul,
    reconstruction_weights,
    verify_mod6_derivative_cir,
)
from treeeval_mv.fixed_weight_closure import (
    fixed_weight_crt_barrier,
    verify_fixed_weight_two_level_closure,
)
from treeeval_mv.factored_vector import (
    FactoredEqualityVectors,
    Overlay,
    apply_overlay,
    query_overlay_direct,
    query_overlay_primewise,
    recursive_type_accounting,
)
from treeeval_mv.gram_search import (
    characteristic_free_size_bound,
    determinant_lift_is_safe,
    equality_product_hpr_family,
    extend_hpr_family_one_label,
    matrix_rank_mod,
    run_canonical_gram_smt,
    verify_gram_characterization,
)
from treeeval_mv.independent_verifier import verify_mv_family_independently
from treeeval_mv.model import (
    MVFamily,
    canonical_customary_set,
    canonical_hpr_offdiagonal_set,
    flip_customary_to_hpr,
)
from treeeval_mv.packed_catalyst import (
    add_vector,
    exhaust_packing_state_map,
    inner_product,
    normalize_arbitrary_catalyst,
    pack_digits,
    packed_bit_length,
    restore_arbitrary_catalyst,
    unpack_digits,
)
from treeeval_mv.projected_selector import (
    has_cross_characteristic_fixed_weight,
    verify_projected_selector_identity,
)
from treeeval_mv.projected_closure_accounting import projected_closure_accounting
from treeeval_mv.succinct_profile_barrier import (
    dynamic_profile_barrier,
    zero_base_scan_report,
)
from treeeval_mv.structured_kernels import (
    factor_matrix_mod,
    kernel_pair_report,
    kernel_pair_to_hpr_family,
    symplectic_oriented_kernels,
)
from treeeval_mv.selector_hitting_set import (
    minimum_selector_hitting_set,
    selector_support,
)
from treeeval_mv.verifier import (
    evaluate_multilinear_polynomial,
    verify_mv_family,
    verify_polynomial_matching_family,
)


class MatchingVectorTests(unittest.TestCase):
    def test_factored_vector_overlay_and_exact_reversal(self) -> None:
        vectors = FactoredEqualityVectors(3, (3, 5, 7))
        explicit = coordinate_equality_hpr_family(
            alphabet_size=3, primes=(3, 5, 7)
        )
        for label in range(vectors.capacity):
            self.assertEqual(vectors.materialize(label), explicit.u[label])
        base = (17, 83, 41)
        for left in range(vectors.capacity):
            for right in range(vectors.capacity):
                overlay = Overlay(left, 19)
                direct = query_overlay_direct(vectors, base, overlay, right)
                independent = query_overlay_primewise(vectors, base, overlay, right)
                self.assertEqual(direct, independent)
            changed = apply_overlay(vectors, base, Overlay(left, 19))
            restored = apply_overlay(vectors, changed, Overlay(left, -19))
            self.assertEqual(restored, base)

    def test_factored_vector_recursive_gate(self) -> None:
        accounting = recursive_type_accounting(
            height=8, label_bits=8, prime_count=8, calls_per_level=4**8
        )
        self.assertEqual(accounting.oracle_leaf_count, 4 ** (8 * 8))
        self.assertEqual(accounting.branch_exponent, 128)
        self.assertGreater(accounting.full_overlay_stack_bits, 8 + 8)

    def test_structured_characteristic_kernel_factorization(self) -> None:
        first, second = symplectic_oriented_kernels(1)
        report = kernel_pair_report(first, second, (3, 5))
        self.assertTrue(report.binary_diagonal_one)
        self.assertTrue(report.hadamard_identity)
        self.assertEqual(report.ranks, (4, 6))
        family = kernel_pair_to_hpr_family(first, second, (3, 5))
        self.assertEqual((family.size, family.dimension), (9, 6))
        self.assertTrue(verify_mv_family(family).ok)
        self.assertTrue(verify_mv_family_independently(family).ok)
        for matrix, prime in ((first, 3), (second, 5)):
            left, right = factor_matrix_mod(matrix, prime)
            for i in range(len(matrix)):
                for j in range(len(matrix)):
                    self.assertEqual(
                        sum(a * b for a, b in zip(left[i], right[j])) % prime,
                        matrix[i][j],
                    )

    def test_coordinate_equality_growing_modulus_family(self) -> None:
        family = coordinate_equality_hpr_family(
            alphabet_size=3, primes=(3, 5, 7)
        )
        self.assertEqual(family.size, 27)
        self.assertEqual(family.dimension, 3)
        self.assertEqual(family.modulus, 105)
        self.assertTrue(verify_mv_family(family).ok)
        self.assertTrue(verify_mv_family_independently(family).ok)
        selector_calibration = coordinate_equality_hpr_family(
            alphabet_size=2, primes=(3, 5)
        )
        self.assertTrue(verify_hpr_selector_identity(selector_calibration).ok)

    def test_branch_aware_torus_barrier_and_projective_calibration(self) -> None:
        barrier = branch_aware_exact_barrier((3, 5))
        self.assertTrue(barrier.exact_selector_impossible)
        self.assertTrue(barrier.nonzero_in_modulus)
        self.assertTrue(barrier.arbitrary_finite_shift_sets)
        self.assertEqual(barrier.telescoping_prime, 3)
        for cycle_length in (3, 5, 7):
            potential = unit_edge_potential(
                cycle_length=cycle_length, modulus=15
            )
            self.assertEqual(len(potential), cycle_length)
            self.assertTrue(verify_unit_edge_potential(potential, 15))
        self.assertEqual(
            multiplicative_unit_subgroups(15),
            (
                (1,),
                (1, 4),
                (1, 11),
                (1, 14),
                (1, 2, 4, 8),
                (1, 4, 7, 13),
                (1, 4, 11, 14),
                (1, 2, 4, 7, 8, 11, 13, 14),
            ),
        )

    def test_fixed_weight_two_level_closure_and_crt_barrier(self) -> None:
        # In one characteristic, phi(x)=x is the mixed-derivative selector.
        positive = verify_fixed_weight_two_level_closure(
            tuple(range(3)), primes=(3,)
        )
        self.assertTrue(positive.ok, positive.errors)
        self.assertEqual(positive.one_level_checks, 36)
        self.assertEqual(positive.two_level_checks, 36**2)
        self.assertEqual(positive.leaf_calls_per_two_levels, 16)

        barrier = fixed_weight_crt_barrier((3, 5))
        self.assertTrue(barrier.impossible)
        self.assertEqual(barrier.modulus, 15)
        self.assertEqual({barrier.idempotent, barrier.complement}, {6, 10})
        self.assertIn("contradicting", barrier.contradiction)

    def test_canonical_sets_mod_6(self) -> None:
        self.assertEqual(canonical_customary_set((2, 3)), frozenset({1, 3, 4}))
        self.assertEqual(canonical_hpr_offdiagonal_set((2, 3)), frozenset({0, 3, 4}))

    def test_dgy_grolmusz_mod_6_calibration(self) -> None:
        polynomial = dgy_grolmusz_polynomial_family(
            primes=(2, 3), exponents=(1, 1), weight=2, universe_size=4
        )
        family = dgy_grolmusz_family(
            primes=(2, 3), exponents=(1, 1), weight=2, universe_size=4
        )
        self.assertEqual(polynomial.size, 6)
        self.assertEqual(polynomial.dimension, 15)
        self.assertLessEqual(
            polynomial.metadata["realized_polynomial_degree"],
            polynomial.metadata["theoretical_degree_bound"],
        )
        self.assertTrue(verify_polynomial_matching_family(polynomial).ok)
        self.assertTrue(verify_mv_family(family).ok)
        self.assertTrue(verify_mv_family_independently(family).ok)

    def test_weight_test_degree_and_truth_table(self) -> None:
        polynomial = weight_congruence_test(
            variable_count=6, prime=2, exponent=2, target_weight=3
        )
        self.assertLessEqual(max(map(len, polynomial), default=0), 3)
        for mask in range(1 << 6):
            point = tuple((mask >> index) & 1 for index in range(6))
            value = evaluate_multilinear_polynomial(polynomial, point, 2)
            expected = 0 if sum(point) % 4 == 3 else 1
            self.assertEqual(value, expected)

    def test_corruption_is_detected_twice(self) -> None:
        family = dgy_grolmusz_family(
            primes=(2, 3), exponents=(1, 1), weight=2, universe_size=4
        )
        corrupted_u = list(family.u)
        corrupted_u[0] = tuple([1] + list(corrupted_u[0][1:]))
        corrupted = MVFamily.build(
            modulus=family.modulus,
            u=corrupted_u,
            v=family.v,
            allowed_offdiagonal=family.allowed_offdiagonal,
            prime_factors=family.prime_factors,
        )
        self.assertFalse(verify_mv_family(corrupted).ok)
        self.assertFalse(verify_mv_family_independently(corrupted).ok)

    def test_hpr_selector_on_odd_modulus_companion(self) -> None:
        customary = dgy_grolmusz_family(
            primes=(3, 5), exponents=(1, 1), weight=2, universe_size=4
        )
        hpr = flip_customary_to_hpr(customary)
        self.assertTrue(verify_mv_family(hpr).ok)
        self.assertTrue(verify_mv_family_independently(hpr).ok)
        report = verify_hpr_selector_identity(hpr)
        self.assertTrue(report.ok, report.errors)

    def test_hpr_selector_rejects_mod_6(self) -> None:
        family = dgy_grolmusz_family(
            primes=(2, 3), exponents=(1, 1), weight=2, universe_size=4
        )
        report = verify_hpr_selector_identity(flip_customary_to_hpr(family))
        self.assertFalse(report.ok)
        self.assertIn("odd", report.errors[0])

    def test_projected_selector_on_odd_modulus_companion(self) -> None:
        customary = dgy_grolmusz_family(
            primes=(3, 5), exponents=(1, 1), weight=2, universe_size=4
        )
        report = verify_projected_selector_identity(flip_customary_to_hpr(customary))
        self.assertTrue(report.ok, report.errors)
        self.assertEqual(report.brute_force_checks, 7200)
        self.assertEqual(report.factor_checks, 7200)

    def test_projected_selector_closure_accounting(self) -> None:
        report = projected_closure_accounting(
            primes=(3, 5), depth=7, label_bits=11
        )
        self.assertEqual(report.modulus, 15)
        self.assertEqual(report.hpr_modulus_bits_per_level, 4)
        self.assertEqual(report.projected_complement_bits, (3, 2))
        self.assertEqual(report.sequential_worst_bits_per_level, 3)
        self.assertEqual(report.hpr_stack_bits, 28)
        self.assertEqual(report.sequential_projected_stack_bits, 21)
        self.assertEqual(report.candidate_outer_oracle_log2_lower_bound, 154)

    def test_selector_exponent_hitting_sets(self) -> None:
        self.assertEqual(selector_support(5, 0, 0), frozenset({0, 1}))
        expected = {3: 2, 5: 3, 7: 4}
        for prime, minimum in expected.items():
            report = minimum_selector_hitting_set(prime)
            self.assertEqual(report.minimum_size, minimum)
            self.assertGreaterEqual(report.minimum_size, report.counting_lower_bound)

    def test_cross_characteristic_fixed_weight_is_impossible(self) -> None:
        self.assertFalse(has_cross_characteristic_fixed_weight(3, 5))
        self.assertFalse(has_cross_characteristic_fixed_weight(5, 3))

    def test_cross_characteristic_holonomy_state_lower_bound(self) -> None:
        forward = verify_holonomy_bound(source_modulus=3, target_modulus=5)
        reverse = verify_holonomy_bound(source_modulus=5, target_modulus=3)
        self.assertEqual(forward.minimum_states, 5)
        self.assertEqual(forward.minimum_bits, 3)
        self.assertEqual(forward.exhaustive_state_counts[:-1], ((1, False), (2, False), (3, False), (4, False)))
        self.assertEqual(forward.exhaustive_state_counts[-1], (5, True))
        self.assertTrue(forward.construction_restores)
        self.assertEqual(reverse.minimum_states, 3)
        self.assertEqual(reverse.exhaustive_state_counts[-1], (3, True))
        self.assertEqual(
            minimum_reversible_holonomy_states(
                source_modulus=6, target_modulus=15, cell_increment=1
            ),
            5,
        )

    def test_exact_gram_rank_and_characteristic_free_bound(self) -> None:
        self.assertEqual(matrix_rank_mod(((1, 1), (1, 1)), 3), 1)
        self.assertEqual(matrix_rank_mod(((1, 1), (1, 2)), 3), 2)
        self.assertTrue(determinant_lift_is_safe(prime=3, rank_bound=2))
        self.assertFalse(determinant_lift_is_safe(prime=3, rank_bound=3))
        self.assertEqual(
            characteristic_free_size_bound(primes=(3, 5), rank_bound=2), 4
        )
        self.assertIsNone(
            characteristic_free_size_bound(primes=(3, 5), rank_bound=3)
        )

    def test_equality_product_gram_family_and_extension(self) -> None:
        family = equality_product_hpr_family(primes=(3, 5), alphabet_size=3)
        self.assertEqual(family.size, 9)
        self.assertEqual(family.dimension, 3)
        report = verify_gram_characterization(family)
        self.assertTrue(report.ok, report.errors)
        self.assertEqual(report.component_ranks, (3, 3))
        self.assertTrue(verify_mv_family(family).ok)
        self.assertTrue(verify_mv_family_independently(family).ok)
        restricted = MVFamily.build(
            modulus=family.modulus,
            u=family.u[:-1],
            v=family.v[:-1],
            allowed_offdiagonal=family.allowed_offdiagonal,
            diagonal_value=1,
            prime_factors=family.prime_factors,
            convention="hpr",
        )
        extension = extend_hpr_family_one_label(restricted)
        self.assertTrue(extension.extended)
        self.assertEqual(extension.family.size, 9)
        self.assertFalse(extend_hpr_family_one_label(family).extended)

    def test_dynamic_profile_state_lower_bound(self) -> None:
        family = equality_product_hpr_family(primes=(3, 5), alphabet_size=2)
        report = dynamic_profile_barrier(family)
        self.assertEqual(report.component_ranks, (2, 2))
        self.assertEqual(report.component_profile_counts, (9, 25))
        self.assertEqual(report.total_profile_count, 225)
        self.assertEqual(report.minimum_bits, 8)
        self.assertEqual(report.exhaustive_components, (True, True))
        scan = zero_base_scan_report(family)
        self.assertEqual(scan.active_candidate_pairs, 16)
        self.assertEqual(scan.update_scalar, 1)
        self.assertTrue(scan.scalar_is_unit)

        hpr = flip_customary_to_hpr(
            dgy_grolmusz_family(
                primes=(3, 5), exponents=(1, 1), weight=2, universe_size=4
            )
        )
        hpr_report = dynamic_profile_barrier(hpr)
        self.assertEqual(hpr_report.component_ranks, (6, 6))
        self.assertEqual(hpr_report.total_profile_count, 15**6)
        self.assertEqual(hpr_report.minimum_bits, 24)

    def test_exact_gram_smt_threshold_mod_15_dimension_2(self) -> None:
        witness = run_canonical_gram_smt(
            primes=(3, 5), dimension=2, size=4, timeout_seconds=30
        )
        no_go = run_canonical_gram_smt(
            primes=(3, 5), dimension=2, size=5, timeout_seconds=30
        )
        self.assertEqual(witness.status, "sat")
        self.assertTrue(witness.verification.ok)
        self.assertTrue(verify_mv_family(witness.family).ok)
        self.assertTrue(verify_mv_family_independently(witness.family).ok)
        self.assertEqual(no_go.status, "unsat")

    def test_exact_gram_smt_normalized_basis_case(self) -> None:
        witness = run_canonical_gram_smt(
            primes=(3, 5),
            dimension=2,
            size=4,
            timeout_seconds=30,
            basis_component=0,
        )
        self.assertEqual(witness.status, "sat")
        self.assertTrue(witness.verification.ok)

        joint = run_canonical_gram_smt(
            primes=(3, 5),
            dimension=3,
            size=9,
            timeout_seconds=30,
            basis_component=0,
            other_basis_labels=(0, 1, 2),
        )
        self.assertEqual(joint.status, "sat")
        self.assertEqual(joint.verification.component_ranks, (3, 3))

    def test_exact_small_maxima_mod_6(self) -> None:
        allowed = canonical_customary_set((2, 3))
        family1, result1, _ = exact_mv_search(
            modulus=6,
            dimension=1,
            allowed_offdiagonal=allowed,
            prime_factors=(2, 3),
        )
        family2, result2, _ = exact_mv_search(
            modulus=6,
            dimension=2,
            allowed_offdiagonal=allowed,
            prime_factors=(2, 3),
        )
        self.assertEqual(result1.maximum_size, 2)
        self.assertEqual(result2.maximum_size, 6)
        self.assertTrue(verify_mv_family(family1).ok)
        self.assertTrue(verify_mv_family(family2).ok)

    def test_structured_constant_weight_maximum(self) -> None:
        sets, result = structured_set_system_search(
            modulus=6,
            universe_size=8,
            weight=6,
            allowed_intersections=canonical_customary_set((2, 3)),
        )
        self.assertEqual(result.maximum_size, 4)
        self.assertEqual(len(sets), 4)

    def test_joint_packed_catalyst_exhaustion(self) -> None:
        report = exhaust_packing_state_map(modulus=3, coordinate_count=4)
        self.assertTrue(report.ok, report.errors)
        self.assertEqual(report.initial_states, 128)
        self.assertEqual(packed_bit_length(3, 4), 7)

    def test_bit_level_packed_catalyst_exhaustion(self) -> None:
        report = exhaust_bit_tape_transitions(modulus=3, coordinate_count=4)
        self.assertTrue(report.ok, report.errors)
        self.assertEqual(report.initial_states, 128)
        self.assertEqual(report.transitions, 1536)

    def test_bit_level_power_of_two_state_count(self) -> None:
        report = exhaust_bit_tape_transitions(modulus=4, coordinate_count=2)
        self.assertTrue(report.ok, report.errors)
        self.assertEqual(report.initial_states, 16)

    def test_packed_updates_and_inner_product(self) -> None:
        digits = (1, 4, 0, 3)
        packed = pack_digits(digits, 5)
        self.assertEqual(unpack_digits(packed, 5, 4), digits)
        updated = add_vector(packed, (4, 2, 1, 3), modulus=5, coordinate_count=4)
        self.assertEqual(unpack_digits(updated, 5, 4), (0, 1, 1, 1))
        self.assertEqual(
            inner_product(updated, (2, 3, 4, 1), modulus=5, coordinate_count=4),
            (0 * 2 + 1 * 3 + 1 * 4 + 1) % 5,
        )

    def test_arbitrary_tape_round_trip(self) -> None:
        for initial in range(128):
            normalized = normalize_arbitrary_catalyst(
                initial, modulus=3, coordinate_count=4
            )
            self.assertIn(normalized.quotient_bit, (0, 1))
            self.assertEqual(restore_arbitrary_catalyst(normalized), initial)

    def test_gf4_and_multiplicity_reconstruction(self) -> None:
        self.assertEqual(f4_mul(2, 2), 3)
        self.assertEqual(f4_mul(2, 3), 1)
        self.assertEqual(len(reconstruction_weights()), 4)
        self.assertTrue(any(weight != F4_ONE for weight in reconstruction_weights()))

    def test_mod6_derivative_cir_kernel(self) -> None:
        family = dgy_grolmusz_family(
            primes=(2, 3), exponents=(1, 1), weight=2, universe_size=4
        )
        report = verify_mod6_derivative_cir(family)
        self.assertTrue(report.ok, report.errors)
        self.assertEqual(report.single_kernel_checks, 72)
        self.assertEqual(report.derivative_mask_checks, 288)
        self.assertEqual(report.pair_kernel_checks, 2592)


if __name__ == "__main__":
    unittest.main()
