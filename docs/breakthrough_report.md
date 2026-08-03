# Initial breakthrough audit: TreeEval via matching vectors

Date: 2026-08-01.

## Bottom line

This project did **not** prove `TreeEval in L`. It finalized the initial
objective in the barrier/pivot sense required by the research prompt:

1. the fixed-modulus explicit-vector target needed for L is unconditionally
   impossible;
2. the exact verifier/search stack is implemented and calibrated twice;
3. the arbitrary-bit catalyst can be packed at the information-theoretic
   `O(d log m)` size, concretely realizing a compression target stated after
   BCKLS Lemma 15;
4. the July sparse-decoding PIR does not meet HPR's recursive CIR interface,
   and the exact failed conversion is identified; and
5. the remaining construction target is narrowed to a batched operational
   representation, with a separate growing-modulus frontier that requires a
   new reconstruction-state mechanism as well as a new MV family.

The construction-first continuation supplies a simple exact growing-modulus
baseline: the coordinate-equality CRT family has `N=a^t`, `d=a`, is directly
logspace uniform, and satisfies HPR's canonical convention.  With the first
`t` odd primes it reaches packed size `Theta(ell log ell)` at constant
alphabet, while exposing HPR's `2^t` oracle branching as the other loss.  It is
not the final object, but it couples the two remaining improvements in a form
that can be generated and verified exactly.

A follow-up exact identity narrows Item 5 further: one CRT output component
can discard its own coefficient-extraction state by using a mixed finite
difference, leaving `O(log(m/p_k)+t)` suspended bits. The obstruction is now
localized to recursively composing different characteristics. A reversible
finite-state repair of that local cross-characteristic obstruction needs at
least a full target field of states, so the surviving transducer must be
nonlocal or change the recursion interface.

The post-email audit adds two sharper barriers. A generic dynamic state that
supports arbitrary MV updates and all Gram queries needs
`prod_i p_i^rank(G_i)` distinguishable profiles, even if its encoding is
nonlinear. Separately, any fixed exponent menu hitting every HPR selector
polynomial has `Omega(p)` entries, so choosing a coefficient from a smarter
menu still retains `Omega(log p)` state. These results rule out two tempting
forms of “succinctness” without claiming a lower bound for trace-aware
recomputation.

## A. Closed route: fixed modulus and explicit coordinates

For fixed `m`, Bhowmick-Dvir-Lovett Theorem 2 plus the now-proved
bounded-torsion PFR theorem gives

`MV(m,d) <= exp(c_m d/log d)`.

HPR requires `N=2^ell`, so `d=Omega_m(ell log ell)`. A classical clean
simulation of the three explicit vector registers therefore takes
`Omega(ell log ell)` bits. Uniformity cannot change the bound. The search for
fixed-`m`, `d=O(ell)` families is retired.

The exact finite no-go ledger is retained only for calibration:

- customary canonical `Z_6`, `d=1`: maximum `N=2`, independently SMT-UNSAT at
  `N>=3`;
- `Z_6`, `d=2`: maximum `N=6`, independently SMT-UNSAT at `N>=7`;
- 6-uniform restricted-intersection families on ground sizes 6, 7, 8, 9:
  maxima 1, 1, 4, 12.

## B. Concrete secondary result: dense catalytic packing

Jointly rank-encode all `D=Theta(d)` ring coordinates as one integer below
`M=m^D`. An arbitrary `ceil(log_2 M)`-bit catalyst is `qM+x` with one quotient
bit `q`. Retain `q`, normalize to `x`, update base-`m` digits by generated long
constants, and restore `qM+x` at the end.

This replaces HPR Remark 2.3's `O(d log(dm))` binary catalyst accounting with
`O(d log m)` and polynomial logspace arithmetic. Two implementations exhaust
all 128 initial tapes for `m=3,D=4`, all coordinates and deltas: 1,536 inverse
transitions each. This improves HPR's exact catalytic parameter but not the
coarse `exp(O(ell^epsilon))` corollary and does not cross the fixed-modulus
dimension barrier.

BCKLS anticipated linear-size high-bit compression after Lemma 15; the result
here is a concrete specialization, not a claim that the target is new.

## C. July PIR: algebra works, recursive CIR does not

The GKS/Gupte query is additive after writing roots in exponent form, and a
pair query tensors to `T^2` server pairs. Exact `GF(4)/Z_6` calibration verifies
all single and pair linear kernel coefficients, including arbitrary derivative
direction masks.

The end-to-end recursion nevertheless fails. GKS answers live in a field of
characteristic `q`, while the query update is in an exponent group of order
coprime to `q`. Coordinatewise one-hot conversion is correct but repeats child
oracle work `d` times per level, yielding `T(h)>=dT(h-1)`. Batching into an
arbitrary-offset field accumulator would require a nontrivial additive map
from a `q`-group to the coprime exponent group; no such homomorphism exists.
Therefore the July server reduction supplies no HPR TreeEval improvement by
itself.

## D. Remaining inside-the-bounds target

A fixed number of *growing* prime factors escapes the fixed-modulus theorem.
BDL's general restricted-family bound permits the tight budget
`d log m=Theta(ell)`. In HPR's convention, existence is exactly equivalent to
zero-one matrices `G_i` with `rank_{F_{p_i}}(G_i)<=d` and
`G_1 circ ... circ G_t=I_N`. If the zero-one determinant lift is safe at every
prime, these ranks lift to rational rank and force `N<=d^t`. Thus the first
unexcluded balanced regime lies near
`d=Omega(sqrt(ell/log ell))` and `log m=O(sqrt(ell log ell))`.

The exact Gram encoding proves `N_max=4` at `d=2,t=2` and `N_max=8` at
`d=2,t=3`, with matching coordinate-equality product constructions. The first
odd-characteristic exception is `d=3,p=3`; the first capacity-improving
instance `(p_1,p_2,d,N)=(3,5,3,10)` remains unresolved after a logged solver
timeout and is not claimed as a no-go.

The coordinate-equality construction extends uniformly to arbitrary `t`:
labels are words in `[a]^t`, component `i` records equality of the `i`th
symbols, and one-hot factorization gives `N=a^t,d=a`.  Taking any prefix gives
exactly `N=2^ell`.  Its packed budget is
`Theta(2^(ell/t) t log t)`, minimized at `Theta(ell log ell)` rather than the
required `Theta(ell)`.  This is now the executable baseline for attempts to
replace equality components by higher-capacity characteristic-dependent
kernels.

Even a family there does not plug into HPR as written: their suspended
reconstruction state costs `h log m`. The naive recomputation shortcut fails
because computing each inner-product state itself must preserve a `log m`-bit
baseline across child calls. The honest growing-modulus target is therefore a
pair of objects:

1. a uniform, near-information-tight canonical MV family; and
2. a transparent, polynomial-time reconstruction-state stack using
   `O(ell+h+log m)` rather than `h log m` bits.

## E. Exact next theorem to seek

The shortest formulation is a batched Gram transducer. It must accept an
implicitly evaluated child value, stream every Gram query used in one HPR
truth-table pass, support balanced updates and exact reversal, and use
`O(ell+h)` live bits in total polynomial time. Static coordinate generation or
a short execution-prefix descriptor is insufficient; both were explicitly
shown to incur a multiplicative `d` or `2^(2ell)` replay at every level.

The projected-selector lemma in `projected_selector_lemma.md` supplies a more
algebraic equivalent target: combine all characteristic-local projected
updates into a recursively closed full-vector update without storing the
complementary CRT residues at every level. A state-independent fixed weight
cannot do this across distinct prime characteristics, by a telescoping
identity proved there. A stateful local repair also fails cheaply: one source
cycle has target holonomy `p in F_q`, so a reversible return-map cycle must
have length divisible by `q`. It therefore needs at least `q` states, restoring
`ceil(log_2 q)` persistent bits. The next theorem must exploit nonlocal batching
across candidates/components or a different recursive interface.

The plug-compatible object is now specified in
`projected_closure_spec.md`: it must return a full same-type `Z_m^d` update,
retain state independent of `log m` through child calls, make only `2^O(t)`
child calls per level, and restore every catalyst exactly. Sequential projected
updates retain `h log(m/p_min)` worst-path state, while candidate-outer
composition has at least `2^(2 ell h)` oracle leaves.

## Verification status

- 28 automated tests pass.
- Primary and independent CRT MV verifiers agree.
- DGY polynomial representation and HPR selector calibration pass where their
  hypotheses apply.
- Z3 independently certifies each logged unstructured threshold no-go.
- All arithmetic is exact; no floating-point computation is used.
- The projected selector passes 7,200 direct global identities and 7,200
  independent factor identities; finite function search confirms the
  cross-characteristic obstruction in both directions between `F_3` and
  `F_5`.
- Exact bit-vector Gram search gives an independently hashed `N=4` SAT / `N=5`
  UNSAT threshold at `(3,5),d=2`; every witness is rechecked by three verifiers.
- Arithmetic and exhaustive permutation checks agree that local reversible
  cross-characteristic holonomy needs 5 states for `F_3 -> F_5` and 3 for
  `F_5 -> F_3`.
- Withdrawn ECCC TR26-044 is excluded.

The author-facing two-page note was reviewed and sent to the HPR authors on
1 August 2026. The packing result remains pending author feedback.
