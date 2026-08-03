# Search retarget after the fixed-modulus barrier

Status: 2026-08-01.

## Retired objective

Do not spend search time looking for a fixed-modulus explicit matching-vector
family with `N=2^ell` and `d=O(ell)`. Bhowmick-Dvir-Lovett Theorem 2, now
unconditional via bounded-torsion PFR, forces
`d=Omega_m(ell log ell)` for every fixed `m`, even without a restricted
off-diagonal avoiding set. SAT/ILP cannot overturn this asymptotic theorem.

The existing `(m,d)=(6,1),(6,2)` maximum-clique exhaustions and constant-weight
set-system exhaustions remain useful as:

- independent calibration of the exact modular encodings;
- bounded no-go data for specific construction classes; and
- tiny test instances for proposed succinct operational representations.

They are not a pipeline whose larger output is expected to reach TreeEval in
`L`.

## New primary object: an operational MV representation

A candidate representation is not scored by explicit vector dimension alone.
For HPR it must supply, exactly and uniformly:

1. a descriptor for every `u_a` and `v_a`, with labels `a in {0,1}^ell`;
2. reversible updates `X <- X + gamma u_a` and
   `Y <- Y + gamma v_a` on a logical evolving state;
3. queries `<X,v_r>` and `<u_r,Y>` modulo `m`;
4. exact inverse operations and restoration of any storage used; and
5. `O(ell+h log m)` live clean bits and polynomial total time.

For a classical `L` algorithm the initial logical vector registers are zero,
so the arbitrary-catalyst condition may be dropped. For a genuine sharpening
of HPR's catalytic theorem it must be retained.

The critical distinction is between a short **static descriptor** for one MV
vector and a short **dynamic descriptor** for a sum of vectors accumulated
along the recursive execution. Logspace coordinate generation solves only the
static problem.

## Exact growing-modulus search object

In HPR's diagonal-one convention, let `G_i` be the component Gram matrix over
`F_{p_i}`. An HPR family exists exactly when

1. every `G_i` is a zero-one matrix with diagonal one;
2. `rank_{F_{p_i}}(G_i) <= d`; and
3. `G_1 circ ... circ G_t = I_N` entrywise.

Necessity is immediate from `G_i=U_i V_i^T`. Conversely, factor every `G_i`
over its field, pad to dimension `d`, and CRT-combine corresponding factor
coordinates. This is a necessary-and-sufficient finite search object, not a
relaxation.

`gram_search.py` encodes the factors in exact bit-vector arithmetic, with
enough width to prevent overflow before modular reduction. Solver witnesses
are reconstructed over `Z_m` and checked by the primary MV verifier, the
independent CRT verifier, and a separate rank/Hadamard verifier.

The first exact thresholds are:

- `t=2,d=2`: `N_max=4` for all odd primes;
- `t=3,d=2`: `N_max=8` for all odd primes.

The lower bound is the coordinate-equality product family `N=d^t`. For `d=2`,
the sharp zero-one determinant bound lifts every odd-characteristic rank bound
to rational rank, and rational Hadamard rank gives the matching upper bound.
The first exceptional odd-characteristic case is therefore `d=3,p=3`: a
nonzero four-by-four zero-one determinant can vanish modulo 3.  The bounded
`(p_1,p_2,d,N)=(3,5,3,10)` instance has now been exhausted by a complete
38-case rank/basis split and is UNSAT.  A positive `N=9` control passes the
same encoding, so the maximum is exactly 9.  The determinant lift extends
this maximum to every pair of distinct odd primes; see
`exceptional_gram_n10_no_go.md`.

## Searchable structural scores

Exact finite search is retargeted to representations whose operations can be
proved, not just compressed bit strings. Candidate families should be logged
with the following scores.

| Score | Finite measurement | Why it matters |
|---|---|---|
| Gram oracle | smallest circuit/decision diagram for `(a,r) -> <u_a,v_r>` | HPR inner products against a sum reduce to repeated oracle calls |
| Update history | maximum number of live terms in each logical register under an exact HPR trace | controls a sparse-sum representation |
| Coefficient alphabet | distinct nonzero coefficients reached per trace level | may allow constant-state cancellation |
| CRT separability | descriptor size when each prime component is represented independently | HPR selectors inspect zero/nonzero residues primewise |
| Recompute recurrence | measured/proved time recurrence when a term label is recomputed rather than stored | separates polynomial replay from quasipolynomial replay |
| Reversible state size | bits needed to invert every update without a hidden transcript | rules out lossy compression masquerading as a data structure |
| Characteristic lift | rank gain over `Q` versus each `F_p`, plus cross-field holonomy state | isolates improvement over `N=d^t` and its recursive cost |

## Structured construction sweep

For restricted-intersection/DGY families, materializing the monomial-evaluation
vector is no longer the goal. Work with the underlying set label and test
whether the following can be computed directly from two labels:

- intersection size modulo each prime;
- the canonical zero-pattern of the matching product;
- the selector coefficient used by HPR; and
- the contribution of a sparse sum of active set labels.

The first three are static and are often easy. The last item is the actual
bottleneck. A headline result must include its update/recompute recurrence.

## Claim rule

No succinct family is claimed from a short Gram oracle alone. A TreeEval
consequence requires an end-to-end simulation of the HPR register interface,
with a proved space bound and a polynomial-time recurrence. This is the guard
against repeating the degree-accounting failure behind withdrawn TR26-044.

The local stateful repair of the projected selector is also retired. Summing a
constant mixed difference around `F_p` creates target-field holonomy `p` in
`F_q`. A reversible state return map must have a cycle length divisible by
`q`, hence at least `q` states. See `cross_characteristic_holonomy.md`. A viable
transducer must batch nonlocally across candidates/components or change HPR's
recursion interface.
