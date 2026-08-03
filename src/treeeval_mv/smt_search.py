"""Independent finite SAT/SMT encoding for bounded no-go cross-checks."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Sequence


@dataclass(frozen=True)
class SMTResult:
    status: str
    threshold: int
    sha256: str
    solver: str
    variables: int
    incompatibility_clauses: int


def clique_threshold_smt2(adjacency: Sequence[int], threshold: int) -> str:
    count = len(adjacency)
    lines = ["(set-logic QF_LIA)"]
    lines.extend(f"(declare-fun x_{index} () Bool)" for index in range(count))
    clauses = 0
    for left in range(count):
        for right in range(left):
            if not (adjacency[left] >> right) & 1:
                lines.append(f"(assert (or (not x_{left}) (not x_{right})))")
                clauses += 1
    cardinality = " ".join(f"(ite x_{index} 1 0)" for index in range(count))
    lines.append(f"(assert (>= (+ {cardinality}) {threshold}))")
    lines.append("(check-sat)")
    lines.append("(exit)")
    return "\n".join(lines) + "\n"


def run_clique_threshold_smt(
    adjacency: Sequence[int], threshold: int, *, timeout_seconds: int = 120
) -> SMTResult:
    solver_path = shutil.which("z3")
    if solver_path is None:
        raise RuntimeError("z3 executable not found")
    encoding = clique_threshold_smt2(adjacency, threshold)
    digest = sha256(encoding.encode("utf-8")).hexdigest()
    clauses = sum(
        1
        for left in range(len(adjacency))
        for right in range(left)
        if not (adjacency[left] >> right) & 1
    )
    with tempfile.TemporaryDirectory(prefix="treeeval-mv-smt-") as directory:
        path = Path(directory) / "instance.smt2"
        path.write_text(encoding, encoding="utf-8")
        process = subprocess.run(
            [solver_path, f"-T:{timeout_seconds}", str(path)],
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout_seconds + 10,
        )
    output = process.stdout.strip().splitlines()
    status = output[0] if output else f"error(exit={process.returncode})"
    version = subprocess.run(
        [solver_path, "-version"], check=False, capture_output=True, text=True
    ).stdout.strip()
    return SMTResult(status, threshold, digest, version, len(adjacency), clauses)
