"""Exact GF(4)/Z_6 calibration of a derivative-PIR-to-CIR lift.

The construction is the smallest multiplicity case of Ghasemi--Kopparty--
Sudan: m=3, characteristic q=2, and matching-vector modulus M=mq=6.
It deliberately works at the linear-kernel level.  Verifying every source and
target pair therefore verifies every database over GF(4), not just sampled
database contents.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Sequence

from .model import MVFamily, Vector


# GF(4) = GF(2)[a]/(a^2+a+1), encoded in polynomial basis as 0,1,a,a+1.
F4_ZERO = 0
F4_ONE = 1
F4_OMEGA = 2


def f4_add(left: int, right: int) -> int:
    return left ^ right


def f4_mul(left: int, right: int) -> int:
    raw = 0
    for bit in range(2):
        if (right >> bit) & 1:
            raw ^= left << bit
    # Reduce a^3 and then a^2 using a^2=a+1.
    if raw & 0b1000:
        raw ^= 0b1110
    if raw & 0b0100:
        raw ^= 0b0111
    return raw & 0b11


def f4_pow(value: int, exponent: int) -> int:
    if exponent < 0:
        return f4_pow(f4_inv(value), -exponent)
    result = F4_ONE
    base = value
    power = exponent
    while power:
        if power & 1:
            result = f4_mul(result, base)
        base = f4_mul(base, base)
        power >>= 1
    return result


def f4_inv(value: int) -> int:
    if value == F4_ZERO:
        raise ZeroDivisionError("zero has no inverse in GF(4)")
    return f4_pow(value, 2)


def _f4_sum(values: Sequence[int]) -> int:
    result = F4_ZERO
    for value in values:
        result = f4_add(result, value)
    return result


def _dot_mod(left: Vector, right: Vector, modulus: int) -> int:
    return sum(a * b for a, b in zip(left, right, strict=True)) % modulus


def _root(exponent: int) -> int:
    return f4_pow(F4_OMEGA, exponent % 3)


def _hasse_monomial(exponent: int, order: int, point: int) -> int:
    if order == 0:
        return f4_pow(point, exponent)
    if order == 1:
        return f4_pow(point, exponent - 1) if exponent & 1 else F4_ZERO
    raise ValueError("the calibration supports multiplicity two only")


OBSERVATIONS = ((0, 0), (0, 1), (1, 0), (1, 1))
CANONICAL_SUPPORT_MOD6 = (0, 1, 3, 4)


def _solve_f4(matrix: list[list[int]], target: list[int]) -> tuple[int, ...]:
    size = len(target)
    augmented = [row[:] + [target[index]] for index, row in enumerate(matrix)]
    for column in range(size):
        pivot = next(row for row in range(column, size) if augmented[row][column])
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        inverse = f4_inv(augmented[column][column])
        augmented[column] = [f4_mul(inverse, value) for value in augmented[column]]
        for row in range(size):
            if row == column or not augmented[row][column]:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                f4_add(value, f4_mul(factor, pivot_value))
                for value, pivot_value in zip(
                    augmented[row], augmented[column], strict=True
                )
            ]
    return tuple(row[-1] for row in augmented)


def reconstruction_weights() -> tuple[int, ...]:
    """Weights mapping value/derivative observations to coefficient Z^0."""

    # observations = E * coefficients, so solve E^T lambda = e_0.
    evaluation = [
        [
            _hasse_monomial(exponent, order, _root(point_exponent))
            for exponent in CANONICAL_SUPPORT_MOD6
        ]
        for point_exponent, order in OBSERVATIONS
    ]
    transpose = [list(column) for column in zip(*evaluation, strict=True)]
    return _solve_f4(transpose, [1, 0, 0, 0])


def _weighted_observations(observations: Sequence[int]) -> int:
    weights = reconstruction_weights()
    return _f4_sum(
        [f4_mul(weight, value) for weight, value in zip(weights, observations, strict=True)]
    )


def single_source_direct_kernel(
    family: MVFamily,
    *,
    target: int,
    source: int,
    random_exponents: Vector,
) -> int:
    """Reconstructed coefficient for one database basis vector."""

    exponent = _dot_mod(family.v[target], family.u[source], 6)
    mask = _root(_dot_mod(random_exponents, family.u[source], 3))
    observations = []
    for point_exponent, order in OBSERVATIONS:
        point = _root(point_exponent)
        observations.append(f4_mul(mask, _hasse_monomial(exponent, order, point)))
    reconstructed = _weighted_observations(observations)
    target_mask = _root(_dot_mod(random_exponents, family.u[target], 3))
    return f4_mul(reconstructed, f4_inv(target_mask))


def _masked_direction_server_term(
    family: MVFamily,
    *,
    target: int,
    source: int,
    random_exponents: Vector,
    point_exponent: int,
    direction_mask: Vector,
) -> int:
    """A first-derivative server term contracted with an arbitrary direction."""

    query = tuple(
        (random_value + point_exponent * target_value) % 3
        for random_value, target_value in zip(
            random_exponents, family.v[target], strict=True
        )
    )
    value = _root(_dot_mod(query, family.u[source], 3))
    direction = _dot_mod(
        tuple(coordinate % 2 for coordinate in family.u[source]), direction_mask, 2
    )
    if not direction:
        return F4_ZERO
    return f4_mul(value, f4_inv(_root(point_exponent)))


def single_source_masked_derivative(
    family: MVFamily,
    *,
    target: int,
    source: int,
    random_exponents: Vector,
    point_exponent: int,
    direction_mask: Vector,
) -> int:
    """Extract the hidden derivative direction by cancelling an arbitrary mask."""

    updated = tuple(
        (mask + target_value) % 2
        for mask, target_value in zip(direction_mask, family.v[target], strict=True)
    )
    with_target = _masked_direction_server_term(
        family,
        target=target,
        source=source,
        random_exponents=random_exponents,
        point_exponent=point_exponent,
        direction_mask=updated,
    )
    without_target = _masked_direction_server_term(
        family,
        target=target,
        source=source,
        random_exponents=random_exponents,
        point_exponent=point_exponent,
        direction_mask=direction_mask,
    )
    return f4_add(with_target, without_target)


def single_source_published_derivative(
    family: MVFamily,
    *,
    target: int,
    source: int,
    random_exponents: Vector,
    point_exponent: int,
) -> int:
    exponent = _dot_mod(family.v[target], family.u[source], 6)
    mask = _root(_dot_mod(random_exponents, family.u[source], 3))
    return f4_mul(
        mask,
        _hasse_monomial(exponent, 1, _root(point_exponent)),
    )


def _pair_server_masked_term(
    family: MVFamily,
    *,
    target: tuple[int, int],
    source: tuple[int, int],
    random_exponents: tuple[Vector, Vector],
    point_exponents: tuple[int, int],
    derivative_orders: tuple[int, int],
    direction_masks: tuple[Vector, Vector],
) -> int:
    result = F4_ONE
    for side in range(2):
        query = tuple(
            (random_value + point_exponents[side] * target_value) % 3
            for random_value, target_value in zip(
                random_exponents[side], family.v[target[side]], strict=True
            )
        )
        result = f4_mul(
            result,
            _root(_dot_mod(query, family.u[source[side]], 3)),
        )
        if derivative_orders[side]:
            direction = _dot_mod(
                tuple(value % 2 for value in family.u[source[side]]),
                direction_masks[side],
                2,
            )
            if not direction:
                return F4_ZERO
            result = f4_mul(result, f4_inv(_root(point_exponents[side])))
    return result


def pair_source_masked_observation(
    family: MVFamily,
    *,
    target: tuple[int, int],
    source: tuple[int, int],
    random_exponents: tuple[Vector, Vector],
    point_exponents: tuple[int, int],
    derivative_orders: tuple[int, int],
    direction_masks: tuple[Vector, Vector],
) -> int:
    """Inclusion-exclusion over adversarial direction masks for a pair query."""

    active_sides = [side for side, order in enumerate(derivative_orders) if order]
    result = F4_ZERO
    for update_bits in product((0, 1), repeat=len(active_sides)):
        masks = [list(direction_masks[0]), list(direction_masks[1])]
        for bit, side in zip(update_bits, active_sides, strict=True):
            if bit:
                masks[side] = [
                    (mask + value) % 2
                    for mask, value in zip(
                        masks[side], family.v[target[side]], strict=True
                    )
                ]
        # In characteristic two every inclusion-exclusion sign is +1.
        result = f4_add(
            result,
            _pair_server_masked_term(
                family,
                target=target,
                source=source,
                random_exponents=random_exponents,
                point_exponents=point_exponents,
                derivative_orders=derivative_orders,
                direction_masks=(tuple(masks[0]), tuple(masks[1])),
            ),
        )
    return result


def pair_source_direct_observation(
    family: MVFamily,
    *,
    target: tuple[int, int],
    source: tuple[int, int],
    random_exponents: tuple[Vector, Vector],
    point_exponents: tuple[int, int],
    derivative_orders: tuple[int, int],
) -> int:
    result = F4_ONE
    for side in range(2):
        exponent = _dot_mod(family.v[target[side]], family.u[source[side]], 6)
        mask = _root(
            _dot_mod(random_exponents[side], family.u[source[side]], 3)
        )
        result = f4_mul(result, mask)
        result = f4_mul(
            result,
            _hasse_monomial(
                exponent,
                derivative_orders[side],
                _root(point_exponents[side]),
            ),
        )
    return result


def _pair_reconstruct(
    family: MVFamily,
    *,
    target: tuple[int, int],
    source: tuple[int, int],
    random_exponents: tuple[Vector, Vector],
    direction_masks: tuple[Vector, Vector] | None,
) -> int:
    weights = reconstruction_weights()
    total = F4_ZERO
    for left_index, (left_point, left_order) in enumerate(OBSERVATIONS):
        for right_index, (right_point, right_order) in enumerate(OBSERVATIONS):
            if direction_masks is None:
                observation = pair_source_direct_observation(
                    family,
                    target=target,
                    source=source,
                    random_exponents=random_exponents,
                    point_exponents=(left_point, right_point),
                    derivative_orders=(left_order, right_order),
                )
            else:
                observation = pair_source_masked_observation(
                    family,
                    target=target,
                    source=source,
                    random_exponents=random_exponents,
                    point_exponents=(left_point, right_point),
                    derivative_orders=(left_order, right_order),
                    direction_masks=direction_masks,
                )
            coefficient = f4_mul(weights[left_index], weights[right_index])
            total = f4_add(total, f4_mul(coefficient, observation))
    target_mask = f4_mul(
        _root(_dot_mod(random_exponents[0], family.u[target[0]], 3)),
        _root(_dot_mod(random_exponents[1], family.u[target[1]], 3)),
    )
    return f4_mul(total, f4_inv(target_mask))


@dataclass(frozen=True)
class DerivativeCIRCalibrationReport:
    ok: bool
    single_kernel_checks: int
    derivative_mask_checks: int
    pair_kernel_checks: int
    errors: tuple[str, ...]


def verify_mod6_derivative_cir(family: MVFamily) -> DerivativeCIRCalibrationReport:
    if family.modulus != 6:
        raise ValueError("the calibration requires an MV family over Z_6")
    dimension = family.dimension
    random_patterns = (
        tuple(0 for _ in range(dimension)),
        tuple(index % 3 for index in range(dimension)),
    )
    direction_patterns = (
        tuple(0 for _ in range(dimension)),
        tuple(index % 2 for index in range(dimension)),
    )
    errors: list[str] = []
    single_checks = 0
    derivative_checks = 0
    pair_checks = 0

    for random_exponents in random_patterns:
        for target in range(family.size):
            for source in range(family.size):
                actual = single_source_direct_kernel(
                    family,
                    target=target,
                    source=source,
                    random_exponents=random_exponents,
                )
                expected = F4_ONE if target == source else F4_ZERO
                single_checks += 1
                if actual != expected:
                    errors.append(
                        f"single kernel target={target}, source={source}: {actual}!={expected}"
                    )
                for point_exponent in (0, 1):
                    for direction_mask in direction_patterns:
                        masked = single_source_masked_derivative(
                            family,
                            target=target,
                            source=source,
                            random_exponents=random_exponents,
                            point_exponent=point_exponent,
                            direction_mask=direction_mask,
                        )
                        published = single_source_published_derivative(
                            family,
                            target=target,
                            source=source,
                            random_exponents=random_exponents,
                            point_exponent=point_exponent,
                        )
                        derivative_checks += 1
                        if masked != published:
                            errors.append(
                                "direction cancellation mismatch "
                                f"target={target}, source={source}, point={point_exponent}"
                            )

    # The pair identity is bilinear. Checking every source basis element proves
    # correctness for every database, while two nontrivial masks test catalyst
    # independence.
    for pattern_index in range(2):
        random_exponents = (
            random_patterns[pattern_index],
            random_patterns[1 - pattern_index],
        )
        direction_masks = (
            direction_patterns[pattern_index],
            direction_patterns[1 - pattern_index],
        )
        for target in product(range(family.size), repeat=2):
            for source in product(range(family.size), repeat=2):
                direct = _pair_reconstruct(
                    family,
                    target=target,
                    source=source,
                    random_exponents=random_exponents,
                    direction_masks=None,
                )
                masked = _pair_reconstruct(
                    family,
                    target=target,
                    source=source,
                    random_exponents=random_exponents,
                    direction_masks=direction_masks,
                )
                expected = F4_ONE if target == source else F4_ZERO
                pair_checks += 1
                if direct != expected or masked != direct:
                    errors.append(
                        "pair kernel mismatch "
                        f"target={target}, source={source}, direct={direct}, masked={masked}"
                    )
    return DerivativeCIRCalibrationReport(
        not errors,
        single_checks,
        derivative_checks,
        pair_checks,
        tuple(errors[:20]),
    )
