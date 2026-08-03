"""Exact matching-vector verification and bounded search tools."""

from .constructions import (
    coordinate_equality_hpr_family,
    dgy_grolmusz_family,
    incidence_family,
)
from .hpr_selector import verify_hpr_selector_identity
from .model import (
    MVFamily,
    PolynomialMatchingFamily,
    canonical_customary_set,
    canonical_hpr_offdiagonal_set,
    flip_customary_to_hpr,
)
from .packed_catalyst import normalize_arbitrary_catalyst, packed_bit_length
from .verifier import verify_mv_family, verify_polynomial_matching_family

__all__ = [
    "coordinate_equality_hpr_family",
    "MVFamily",
    "PolynomialMatchingFamily",
    "canonical_customary_set",
    "canonical_hpr_offdiagonal_set",
    "dgy_grolmusz_family",
    "flip_customary_to_hpr",
    "incidence_family",
    "normalize_arbitrary_catalyst",
    "packed_bit_length",
    "verify_hpr_selector_identity",
    "verify_mv_family",
    "verify_polynomial_matching_family",
]
