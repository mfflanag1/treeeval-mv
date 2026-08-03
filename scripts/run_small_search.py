#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
import time

from treeeval_mv.clique_search import exact_mv_search, structured_set_system_search
from treeeval_mv.constructions import incidence_family
from treeeval_mv.independent_verifier import verify_mv_family_independently
from treeeval_mv.model import canonical_customary_set
from treeeval_mv.smt_search import run_clique_threshold_smt
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


def main() -> None:
    allowed = canonical_customary_set((2, 3))
    ledger: dict[str, object] = {
        "research_status": {
            "asymptotic_target": "fixed-modulus explicit d=O(ell)",
            "status": "retired: ruled out by BDL Theorem 2 plus bounded-torsion PFR",
            "current_use": "finite calibration and structural-complexity measurement only",
            "retarget_document": "docs/search_retarget.md",
        },
        "scope": {
            "modulus": 6,
            "prime_factors": [2, 3],
            "offdiagonal_set": sorted(allowed),
            "convention": "customary diagonal-zero canonical-set MV",
        },
        "unstructured_exact": [],
        "structured_constant_weight": [],
    }
    for dimension in (1, 2):
        started = time.monotonic()
        family, clique, adjacency = exact_mv_search(
            modulus=6,
            dimension=dimension,
            allowed_offdiagonal=allowed,
            prime_factors=(2, 3),
        )
        smt = run_clique_threshold_smt(adjacency, clique.maximum_size + 1)
        primary = verify_mv_family(family)
        independent = verify_mv_family_independently(family)
        ledger["unstructured_exact"].append(
            {
                "dimension": dimension,
                "maximum_N": clique.maximum_size,
                "bounded_no_go": f"no family with N >= {clique.maximum_size + 1}",
                "candidate_vertices": clique.vertices,
                "compatibility_edges": clique.edges,
                "branch_and_bound_nodes": clique.search_nodes,
                "wall_seconds": round(time.monotonic() - started, 6),
                "maximum_clique_proof": "exact branch-and-bound with greedy-color upper bounds",
                "independent_smt_threshold_check": {
                    "status": smt.status,
                    "threshold": smt.threshold,
                    "encoding_sha256": smt.sha256,
                    "solver": smt.solver,
                    "variables": smt.variables,
                    "incompatibility_clauses": smt.incompatibility_clauses,
                },
                "witness": family_payload(family),
                "witness_primary_verified": primary.ok,
                "witness_independent_crt_verified": independent.ok,
            }
        )
    for universe_size in range(6, 10):
        started = time.monotonic()
        sets, clique = structured_set_system_search(
            modulus=6,
            universe_size=universe_size,
            weight=6,
            allowed_intersections=allowed,
        )
        family = incidence_family(
            modulus=6,
            sets=sets,
            universe_size=universe_size,
            allowed=allowed,
            primes=(2, 3),
        )
        primary = verify_mv_family(family)
        independent = verify_mv_family_independently(family)
        ledger["structured_constant_weight"].append(
            {
                "universe_size": universe_size,
                "weight": 6,
                "candidate_sets": clique.vertices,
                "maximum_N": clique.maximum_size,
                "bounded_no_go_within_class": f"no 6-uniform incidence family with N >= {clique.maximum_size + 1}",
                "witness_sets_zero_based": [sorted(values) for values in sets],
                "witness_primary_verified": primary.ok,
                "witness_independent_crt_verified": independent.ok,
                "wall_seconds": round(time.monotonic() - started, 6),
            }
        )
    target = ROOT / "results" / "no_go_ledger.json"
    target.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(target)


if __name__ == "__main__":
    main()
