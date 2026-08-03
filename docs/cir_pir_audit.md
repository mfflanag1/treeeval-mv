# July 2026 sparse-decoding PIR versus HPR's CIR interface

Status: exact small-case algebraic calibration passed on 2026-08-01, but the
end-to-end lift fails the polynomial-time recursion audit. No TreeEval
improvement is claimed.

## 1. Direct-interface verdict

The Gupte-Ragavan result supplies optimally sparse decoding polynomials and
then invokes the derivative matching-vector PIR of
Ghasemi-Kopparty-Sudan (GKS). That published PIR is **not directly a CIR
scheme** in the sense of HPR Definition 4.1.

| HPR obligation | Published derivative PIR | Verdict |
|---|---|---|
| Additively masked deterministic query | Write each root query in exponent form: `beta*b_j^v_a` becomes `r+c_j v_a mod m` | yes |
| Query a pair `a||b` from separate access to `a` and `b` | The PIR is stated for one index | tensor two copies; costs `T^2` server pairs |
| Perfect correctness for every catalyst state | GKS correctness is algebraic and probability one | yes after exponent normalization |
| Query and reconstructed answer have the same type | Queries are exponent vectors; answers are field derivative vectors | no |
| Reconstruction state is short and computable through the child oracle | Chain-rule contraction needs the full hidden direction `v_a mod q` | no as stated |
| Server work streams over the database in low space | Every value/directional-derivative sum is a database scan | yes |
| Reconstruction is a sum of independently processed server answers | Multiplicity interpolation is linear over the field | yes for field-valued records |

The two failures are exactly the obstacles HPR identify after their Section
4.4: type mismatch and recursive composability. Merely substituting a
`T`-sparse polynomial for HPR's inclusion-exclusion is therefore invalid.

## 2. Partial white-box lift that survives the small exact test

There is a way to remove the hidden-direction failure without materializing a
direction vector or a tensor product. It does **not** remove the type failure.

Assume the multiplicity-two GKS setting:

- `q` is the characteristic of a finite field `F`;
- `gcd(q,m)=1`, `M=qm`, and `F` contains the `m`th roots of unity;
- `(u_i,v_i)` is the required matching-vector family over `Z_M^d`;
- `B={b_1,...,b_T}` is the multiplicity-two interpolation set obtained from a
  `T`-sparse decoding polynomial.

Represent a query point by an exponent vector in `Z_m^d`. For target `a` and
server `j`, the deterministic query update is

`R <- R + c_j (v_a mod m)`, where `b_j=omega^c_j`.

The derivative direction is `v_a mod q`. Keep a separate catalytic direction
register `W in F_q^d`, whose initial contents may be adversarial. For a streamed
server sum `H(W)` that is linear in this direction,

`H(W+v_a)-H(W)`

is exactly the required directional derivative and is independent of the
initial `W`. For a pair query and the mixed derivative, use the four-term
identity

`H(W+v_a,W'+v_b)-H(W+v_a,W')-H(W,W'+v_b)+H(W,W')`.

Each term is computed by another streaming pass over the database. The child
oracle performs and reverses the `v_a,v_b` updates, so every catalytic register
is restored. No `d^2` Hessian, direction vector, or transcript is stored.

One can apparently remove the field/exponent type mismatch coordinatewise: to
retrieve a record coordinate whose alphabet is a fixed finite set `A`, run the
linear PIR on the indicator database `[DB_i[c]=rho]` for each `rho in A`.
Perfect correctness reconstructs a one-hot answer, which can be decoded in
constant free space for fixed parameters.

The fatal point is recursion. The direction-mask updates call the child oracle.
If output coordinates are processed one at a time, those child calls are
repeated `d` times, so the running time obeys

`T(h) >= d T(h-1)`.

For every growing matching-vector dimension and `h=Theta(ell)`, this is
superpolynomial. Storing all `d` reconstructed field coordinates in a catalytic
accumulator amortizes the child calls, but leaves an arbitrary additive offset.
Removing that offset while converting the field record into an exponent update
would require a nontrivial additive homomorphism from the characteristic-`q`
answer group to the exponent group of order coprime to `q`; every such
homomorphism is trivial. A nonlinear one-hot decode again needs access to the
before/after accumulator coordinate and restores the `d` replays.

The random root mask is removed in the same style. An arbitrary catalytic
register `U` gives

`<R,U+u_a>-<R,U> = <R,u_a>`.

Thus the child oracle supplies the target-dependent quantity while the
adversarial initial register cancels.

The direction-mask portion is slightly more general than HPR Definition 4.1 because
`AnswerAndReconstruct` makes a constant number of child-oracle calls. Those
calls are balanced and reversible for one field-valued output. They are not
constant in number for an entire exponent-vector output after type conversion.

## 3. Exact calibration

The implementation in `src/treeeval_mv/gks_cir_calibration.py` uses the
smallest nontrivial instance:

- `m=3`, `q=2`, `M=6`, and `F=GF(4)`;
- the independently verified six-element DGY/Grolmusz-style MV family over
  `Z_6`, dimension 15;
- interpolation points `{1,omega}` with value and first Hasse derivative;
- exact reconstruction weights `(omega+1, omega+1, omega, omega+1)` in the
  integer encoding `(3,3,2,3)` of `GF(4)`.

The test checks:

- all 72 single-index source/target kernel coefficients under two exponent
  masks;
- 288 comparisons between the published derivative and arbitrary-direction
  mask cancellation; and
- all 2,592 pair source/target kernel coefficients under two exponent-mask and
  two direction-mask configurations.

The pair check is at the linear-kernel level, so it proves the finite identity
for every `GF(4)` database, not merely sampled database contents. The direct
derivative implementation and the masked inclusion-exclusion implementation
agree exactly.

## 4. Parameter consequence: none without a same-type decoder

For the Gupte-Ragavan optimal decoder with `T=k+1`, GKS uses a matching-vector
modulus with `k+1=T` prime factors and hence has the attractive dimension

`d=exp(O(ell^(1/T) (log ell)^(1-1/T)))`.

The algebraic pair kernel uses only `T^2` server pairs, versus HPR's `2^(2T)`
principal server branches at the same number of prime factors. But the
field-to-exponent conversion changes the actual recursion to at least
`d*T^2` child work per level. Thus the attractive server count does not produce
a polynomial-time TreeEval algorithm.

Gupte-Ragavan's fixed-`T` constructions are unconditional for `T<=15` and
conditional on their number-theoretic conjecture in general. To affect HPR,
one additionally needs a decoder whose queries and linearly reconstructed
answers live in the same reversible update group, or a batched type-conversion
primitive that avoids both an arbitrary-offset accumulator and `d` recursive
replays. The July result supplies neither.

Primary sources:

- [HPR, ECCC TR26-022](https://eccc.weizmann.ac.il/report/2026/022/)
- [GKS, arXiv:2411.11611](https://arxiv.org/abs/2411.11611)
- [Gupte-Ragavan, arXiv:2607.22033](https://arxiv.org/abs/2607.22033)
