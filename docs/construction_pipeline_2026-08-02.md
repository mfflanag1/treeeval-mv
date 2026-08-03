# Construction pipeline: first 1 -> 3 pass

Date: 2026-08-02.

## 1. Characteristic-dependent kernels

The new Gram-level pipeline accepts two binary matrices, computes their exact
ranks in their respective characteristics, factors each matrix, CRT-lifts the
factors to explicit vectors, and submits the result to both MV verifiers.

The first structured candidate uses the alternating form `B` over
`F_3^(2n)`:

- `G[x,y]=1[B(x,y)=0]`;
- `H[x,y]=1[x=y or B(x,y)=1]`.

It is a canonical pair because every off-diagonal one of `G` is a zero of
`H`.  At `n=1`, it gives `N=9` with component ranks `(4,6)` over `(F_3,F_5)`.
The explicit lifted family passes both verifiers.  It does not beat the
coordinate-equality family, which has `N=9,d=3`.  At the next scale, the best
resonant rank found for the second component is 60 at `N=81`, again not
competitive.  The candidate remains in the generator because changing the
oriented relation or the geometry now requires only a new matrix constructor;
factorization and verification are automatic.

## 2. Factored recursive-vector prototype

`FactoredEqualityVectors` stores a label and generates every vector coordinate
directly from its word.  A logical overlay

`base + scalar * u_label`

can be queried without materializing `u_label`.  Two implementations are
checked independently:

1. direct CRT-coordinate generation and modular inner product;
2. primewise equality evaluation followed by CRT.

For the `N=27,d=3,m=105` calibration they agree on all 729 label pairs for an
arbitrary nonzero base.  Applying every tested overlay and then its negative
restores the base exactly.

This is a valid one-level changed vector type.  It is not yet a recursively
closed one: nested calls either retain one descriptor per level or materialize
the overlay and retain the parent's reconstruction state.

## 3. Acceptance gates

The current verdict is deliberately split:

- exact MV algebra: **pass**;
- two independent MV verifiers: **pass**;
- implicit-versus-explicit coordinates: **pass**;
- arbitrary-base query agreement: **pass**;
- exact balanced reversal: **pass**;
- polynomial-time recursive composition for growing `t`: **fail**.

HPR's two-child selector makes `4^t` child-oracle calls per level.  At
`h=t=ell`, the oracle-leaf exponent is `2 ell^2`, whereas the input log-scale
is `Theta(ell)`.  The polynomial degree required therefore grows with `ell`.
Replacing `4^t` by merely `poly(t)` would still give `t^h`, which is
superpolynomial in the balanced regime.  The needed vector type must batch
the child-oracle calls across a level or have constant branching; on-demand
coordinate generation alone is insufficient.

Run `PYTHONPATH=src python3 scripts/run_construction_pipeline.py` to reproduce
the machine-readable report.

## Next construction move

The next kernel sweep should search algebraic pairs where both components are
characteristic-dependent, rather than pairing one anomalous kernel with an
ordinary complement or orientation.  A candidate advances only if

`log N / (d log m)` stays bounded below at growing scales.

For the recursive type, the next positive target is a batched child-update
primitive: compute a child matching vector once into a reversible succinct
handle, apply all scalar multiples required by a parent, and erase the handle
without retaining an `ell`- or `log m`-bit descriptor at every ancestor.
