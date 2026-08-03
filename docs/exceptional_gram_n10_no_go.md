# Exact dimension-3 maximum for two odd CRT primes

Date: 2 August 2026.

## Result

For canonical HPR matching-vector families over `Z_(pq)` in dimension 3,
where `p,q` are any two distinct odd primes, the maximum family size is
exactly 9.  The lower bound is the coordinate-equality product family.  The
only characteristic-dependent upper-bound case reduces to an exhaustive
UNSAT proof for `N=10` over `Z_15`.

This is a bounded no-go, not an asymptotic matching-vector upper bound and not
a TreeEval-in-L algorithm.

## Exhaustive rank split

Let `G_3,G_5` be the binary component Gram matrices.  They have diagonal one,
their off-diagonal supports are disjoint, and their Hadamard product is
`I_N`.  Each has rank at most 3 in its own characteristic.

If both ranks are at most 2, every 3-by-3 binary minor that vanishes modulo 3
or 5 vanishes over the rationals by the exact zero-one determinant bound.
Thus both rational ranks are at most 2, while

`N = rank_Q(I_N) = rank_Q(G_3 hadamard G_5) <= 2*2 = 4`.

Hence an `N=10` witness has rank 3 in at least one component.  This leaves:

- rank `(3,<=2)` and `(<=2,3)`, one normalized factor case each; and
- rank `(3,3)`.  Normalize a `U_3` basis containing label 0 into labels
  `0,1,2`, sort the remaining seven labels by their `U_3` rows, and enumerate
  the other two labels completing `u_{5,0}` to a `U_5` basis.  There are
  `C(9,2)=36` cases.

All 38 cases are UNSAT in the exact finite bit-vector encoding.  The encoding
uses widths large enough to prevent overflow before modular reduction.  Each
case has its own SHA-256 hash in `results/exceptional_gram_n10_bv.json`.

As a positive control, the identical joint-basis encoding at `N=9`, with
both bases on labels `0,1,2`, returns SAT and reconstructs a rank-`(3,3)`
family that passes the Gram verifier.  This guards against an accidentally
overconstrained symmetry split.

An independent `QF_NIA` encoding is also provided by the same generator and
can be logged separately when run.  A spot-check proved the basis-`(0,1,2)`
case UNSAT; the remaining attempted integer cases timed out and are not used
as evidence.  No timeout is interpreted as a no-go.

## Reduction from all odd prime pairs

For a 4-by-4 zero-one matrix, every determinant has absolute value at most 3.
Consequently, for every prime `q>=5`, a binary matrix has rank at most 3 over
`F_q` if and only if it has rank at most 3 over the rationals.  If both CRT
primes are at least 5, the rational Hadamard-rank inequality immediately gives
`N<=3*3=9`.

The only remaining prime pair has the form `(3,q)` with `q>=5`.  The possible
binary rank-at-most-3 matrices in the second component are independent of
`q`, by the preceding determinant observation; they are exactly the matrices
represented by the mod-5 component of the exhaustive search.  The 38-case
`(3,5)` exhaustion therefore rules out `N=10` simultaneously for every such
`q`.  Together with the product lower bound, this proves the stated result.

## Reproduction

Run:

```text
PYTHONPATH=src python scripts/run_exceptional_gram_exhaustion.py --encoding bv
PYTHONPATH=src python scripts/run_exceptional_gram_exhaustion.py --encoding nia
```

The generator writes one ledger per encoding under `results/` and declares
the maximum only when every case returns `unsat`.
