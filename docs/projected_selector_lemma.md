# Characteristic-local projected selector

Status: exact one-level lemma, verified 2026-08-01. This is not by itself a
recursive TreeEval algorithm.

## 1. Statement

Let `m=prod_i p_i` be a product of distinct odd primes and use HPR's
diagonal-one canonical MV convention. Fix a distinguished component `k`.
For a target pair `(a,b)` and a candidate pair `(r,s)`, write

`A_i=<u_a,v_r> mod p_i`, `B_i=<u_b,v_s> mod p_i`.

Every `A_i,B_i` is in `{0,1}`. Both vectors of residues are all ones exactly
at the target; a non-target supplies a zero in at least one component.

For `i != k`, retain HPR's coefficient selector. Given base inner products
`g1_i,g2_i`, choose a nonzero coefficient `alpha_i` and its exponent `beta_i`
in

`sum_{c,d in {0,1}} (-1)^(c+d) X^((g1_i+c)(g2_i+d))`.

For component `k`, do **not** choose a monomial. Under shift bits `c_k,d_k`,
use the field value

`q_k=(g1_k+c_k A_k)(g2_k+d_k B_k) mod p_k`

as the weight. At every other component, include the term only when the
corresponding shifted product equals `beta_i`. Multiply by
`prod_{i != k} alpha_i^(-1)` in `F_{p_k}`, by the usual sign over all shift
bits, and by the CRT idempotent that is one modulo `p_k` and zero at the other
primes.

The total coefficient of `(r,s)`, reduced modulo `p_k`, is

`A_k B_k * prod_{i != k}(alpha_i) * prod_{i != k}(alpha_i)^(-1)`.

It is therefore one at `(a,b)` and zero everywhere else. All other CRT
components receive zero.

The distinguished factor is just the exact mixed finite difference

`xy-(x+1)y-x(y+1)+(x+1)(y+1)=1`.

## 2. State consequence

Run the inner-product routine of HPR Lemma 3.8 modulo `m/p_k`, rather than
modulo `m`. Its three temporary inner products then occupy only
`O(log(m/p_k))` bits across the four oracle calls and still recover every
`g1_i,g2_i` needed for `i != k`.

The one-level algorithm therefore never needs
`g1_k,g2_k,alpha_k,beta_k` across an oracle call. Its suspended reconstruction
state is

`O(t + sum_{i != k} log p_i) = O(t + log(m/p_k))`

rather than `O(log m)`. Arithmetic modulo `p_k` is needed only during the
truth-table scan, when no recursive oracle is active.

This localization is potentially useful for an asymmetric modulus with one
dominant prime. It does **not** immediately improve HPR's recursive theorem:
an oracle query changes a full matching vector in every CRT component, while a
projected child invocation changes only one component. Sequentially invoking
all projections restores the worst complementary state unless a new
cross-component transducer is supplied.

## 3. Why the obvious all-components lift is impossible

One might try to replace coefficient extraction at a source prime `p` by a
fixed weight `phi:F_p -> F_q` while accumulating output in a different
characteristic `q`. Such a weight would have to satisfy, for every `x,y`,

`phi(xy)-phi((x+1)y)-phi(x(y+1))+phi((x+1)(y+1))=1` in `F_q`.

For distinct odd primes `p,q`, this is impossible. Choose
`y notin {0,-1}` and sum over every `x in F_p`. Multiplication by `y` and by
`y+1` are permutations, so all four sums cancel. The left side is zero; the
right side is `p`, which is nonzero in `F_q`.

This is the exact characteristic barrier behind the failed naive product of
lifted residue weights. It also explains why HPR's Boolean coefficient
indicators, whose meanings survive CRT, are necessary for the non-output
components.

## 4. Verification

`src/treeeval_mv/projected_selector.py` contains two checks:

1. a direct global sum over every shift pattern; and
2. an independent factor-by-factor mixed-difference calculation.

On the `Z_15` companion family the verifier exhausts both choices of output
prime, all 225 base-inner-product pairs, and all 16 pairs of canonical delta
patterns: 7,200 global identities and 7,200 factor identities. Exhaustive
function search independently confirms the cross-characteristic no-go for
`F_3 -> F_5` and `F_5 -> F_3`.

## 5. Revised breakthrough target

The missing object can now be stated more sharply than “compress every HPR
state”: construct a recursively closed transducer that turns the collection
of characteristic-local projected updates into a full matching-vector update
without retaining `log(m/p_k)` complementary information per level. A proof
that no such transducer exists in a natural black-box model would also be a
useful barrier theorem.
