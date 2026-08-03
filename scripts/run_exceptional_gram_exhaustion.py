#!/usr/bin/env python3
"""Exhaust the first characteristic-dependent HPR Gram case."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import combinations
import json
from pathlib import Path

from treeeval_mv.gram_search import run_canonical_gram_smt


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--encoding", choices=("bv", "nia"), default="bv")
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    cases: list[tuple[str, dict[str, object]]] = [
        ("rank_3_le_2", {"basis_component": 0, "other_rank_bound": 2}),
        ("rank_le_2_3", {"basis_component": 1, "other_rank_bound": 2}),
    ]
    cases.extend(
        (
            f"rank_3_3_basis_0_{left}_{right}",
            {
                "basis_component": 0,
                "other_basis_labels": (0, left, right),
            },
        )
        for left, right in combinations(range(1, 10), 2)
    )

    def run(case: tuple[str, dict[str, object]]) -> dict[str, object]:
        name, options = case
        result = run_canonical_gram_smt(
            primes=(3, 5),
            dimension=3,
            size=10,
            timeout_seconds=args.timeout,
            encoding=args.encoding,
            solver_seed=0,
            **options,
        )
        return {
            "case": name,
            "options": options,
            "status": result.status,
            "encoding_sha256": result.encoding_sha256,
            "solver": result.solver,
            "wall_seconds": round(result.wall_seconds, 6),
            "variables": result.variables,
            "constraints": result.constraints,
        }

    reports: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(run, case) for case in cases]
        for future in as_completed(futures):
            report = future.result()
            reports.append(report)
            print(f"{report['case']}: {report['status']}", flush=True)
    reports.sort(key=lambda report: str(report["case"]))
    all_unsat = all(report["status"] == "unsat" for report in reports)
    payload = {
        "claim": "maximum canonical HPR family size over Z_15 in dimension 3 is 9",
        "lower_bound": "the coordinate-equality product family has N=9",
        "upper_bound_target": "N=10 is impossible",
        "exhaustive_split": {
            "rank_argument": (
                "if both binary component ranks are at most 2, determinant lift "
                "and the rational Hadamard rank inequality give N<=4"
            ),
            "rank_3_le_2_cases": 2,
            "rank_3_3_cases": 36,
            "rank_3_3_basis_argument": (
                "normalize a U basis containing label 0 in the mod-3 component; "
                "enumerate the two additional basis labels containing label 0 "
                "in the mod-5 component"
            ),
        },
        "encoding": args.encoding,
        "all_cases_unsat": all_unsat,
        "conclusion": "maximum N is exactly 9" if all_unsat else "incomplete",
        "cases": reports,
    }
    target = ROOT / "results" / f"exceptional_gram_n10_{args.encoding}.json"
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(target)


if __name__ == "__main__":
    main()
