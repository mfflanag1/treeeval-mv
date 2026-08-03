# Parameter map and fixed-modulus barrier

Status checked: **2026-08-01**.  Symbols in this note are `N` for MV-family
size, `d` for vector dimension, `m` for modulus, `t` for the number of prime
factors of `m`, and `ell` for HPR's label length, so `N=2^ell`.

## 1. Exact HPR search specification

HPR Theorem 3.1 assumes distinct odd primes `p_1,...,p_t`, their product
`m`, and an `O(ell)`-space-uniform HPR-convention family of size `2^ell` in
`Z_m^d`.  It obtains

- free space `O(ell + h log m)`;
- catalytic space `O(d log(d m))` on an arbitrary bit-initialized catalyst;
- time `poly(2^(ell+h t))`;
- the mild size condition `d log m <= poly(2^(h+ell))`.

HPR's convention has inner products in the `{0,1}` CRT cube, diagonal equal
to one, and off diagonal unequal to one.  HPR Remark 2.11 converts this to the
customary convention used by DGY/BDL by appending one coordinate:
`u'_i=(u_i,1)`, `v'_j=(-v_j,1)`.  The customary diagonal is zero and every
off-diagonal product lies in the nonzero canonical set.

To obtain a **polynomial-time classical-logspace** algorithm merely by storing
HPR's explicit vector registers as ordinary clean workspace, the finite search
target would have to satisfy all of the following:

1. `m` is a fixed product of a fixed number `t` of distinct odd primes.  Fixed
   `t` makes the recursion polynomial time; fixed `log m` is forced by the
   `h log m` free-space term on instances with `h=Theta(ell)`.
2. `N=2^ell`.
3. `d log m=O(ell)`, hence `d=O(ell)` for fixed `m`.
4. The family is logspace uniform in HPR Definition 2.10: on
   `(1^N,i,j,p_1,...,p_t)`, output `u_i,v_j` in
   `O(log N+log d+log m)=O(ell)` space.
5. Exact canonical inner-product conditions hold.

The arbitrary-catalyst representation trick costs `O(d log(d m))`, but a
standard L machine starts with clean workspace and need not pay that validity
offset (compare HPR Remark 2.3 and footnote 6).  Thus `d=O(ell)`, not
`O(ell/log ell)`, is the correct optimistic explicit-vector target.

That target is impossible at fixed modulus under the known upper bound below.
Uniformity cannot help: uniform families are a subclass of all MV families.

## 2. One-table parameter comparison

| Source/object | Modulus and avoiding set | Size versus dimension | Uniformity | Consequence for HPR |
|---|---|---|---|---|
| HPR Thm. 3.1, exact interface | `m=prod_i p_i`, odd distinct primes; canonical CRT pattern after convention flip | Requires `N=2^ell`; accepts arbitrary `d` satisfying the mild polynomial-size condition | `O(ell)`-space uniform | Free `O(ell+h log m)`, catalyst `O(d log(dm))`, time `poly(2^(ell+h t))` |
| Grolmusz / BBR / DGY construction, fixed `t` (HPR Cor. 2.13(1)) | Fixed squarefree composite `m`; canonical set | `N=2^ell`, `d=exp(O(ell^(1/t) (log ell)^(1-1/t)))` | HPR Appendix A proves logspace uniformity | Polynomial time and `O(log n)` free space, but subpolynomial catalyst |
| Exact DGY Lemma 11 form | `m=prod_i p_i`; choose `p_i^e_i>w^(1/t)`, `h_0>=w` | `N=binom(h_0,w)` and vector dimension `sum_{j<=D} binom(h_0,j)`, `D=max_i p_i^e_i` | Coordinate enumeration is logspace uniform (HPR Appendix A) | This is the implemented calibration/generator |
| Efremenko MV/LDC use | Notably `m=511=7*73`; same canonical MV construction plus a 3-sparse decoding polynomial | Does not improve the Grolmusz/DGY MV dimension; improves query/server count | Explicit for its parameters | Useful for branching/PIR refinements, not an escape from the dimension barrier |
| DGY unconditional general upper bound (Thm. 30) | Any `m`, any `S` excluding zero; for every prime `p|m` | `N <= 5 m^d / p^((d-1)/2)` | Applies to all families | For fixed `m`, only forces `d=Omega(ell)`; by itself it leaves the L target open |
| Bhowmick–Dvir–Lovett Thm. 2 + bounded-torsion PFR | Fixed `m`, even with unrestricted nonzero off diagonal | `N <= exp(c_m d/log d)` | Applies to all families | `N=2^ell` forces `d=Omega(ell log ell)`, ruling out the explicit-vector L target |
| Lower-bound-compatible hypothetical in HPR Remark 3.4 | Fixed constant `m,t` | `N=2^ell`, `d=O(ell log ell)` | Must be sufficiently uniform | Gives polynomial time in `O(log n log log n)` classical space; HPR does **not** say this gives L |
| Gupte–Ragavan TR26-126 (2026-07-24) | Special products of primes plus sparse `S`-decoding polynomials | Optimal decoding-polynomial sparsity `t+1` (conditional in general; unconditional for up to 15 servers); the underlying MV dimension remains Grolmusz/DGY | Explicit finite instances; asymptotic result has number-theoretic assumptions | Potentially reduces servers/recursive branches.  It does not reduce `d` and therefore does not remove the storage barrier alone |

Primary links: [HPR / ECCC TR26-022](https://eccc.weizmann.ac.il/report/2026/022/),
[DGY paper](https://www.cs.princeton.edu/~zdvir/papers/DvirGopalanYekhanin10.pdf),
[BDL paper](https://arxiv.org/abs/1204.1367),
[bounded-torsion PFR theorem](https://arxiv.org/abs/2404.02244), and
[Gupte–Ragavan / arXiv:2607.22033](https://arxiv.org/abs/2607.22033).

## 3. Barrier derivation

Let a fixed-modulus customary MV family have `N=2^ell` and dimension `d`.
BDL gives, now unconditionally for fixed bounded-torsion groups,

`2^ell <= exp(c_m d/log d)`.

Taking logarithms gives `ell=O_m(d/log d)`, hence
`d=Omega_m(ell log ell)`.  HPR's three vector registers therefore require
`Omega(ell log ell)` clean bits even before any arbitrary-catalyst encoding
overhead.  Since `log(input length)=Theta(h+ell)` and TreeEval includes
instances with `h=Theta(ell)`, this exceeds logarithmic space by a
`Theta(log ell)` factor.

This is a barrier to the **ordinary explicit fixed-modulus MV instantiation of
HPR Theorem 3.1**.  It is not a lower bound for TreeEval, catalytic information
retrieval, succinct/on-demand representations, or non-MV algorithms.

## 4. What the verifier actually checks

The implementation separates three properties that are easy to conflate:

1. Customary MV: diagonal zero, off diagonal in the prescribed set, by exact
   arithmetic modulo `m`.
2. DGY polynomial matching: `f_i(x_i)=0`, `f_i(x_j)` lies in the canonical
   avoiding set for `i!=j`, and coefficient/evaluation vectors reproduce every
   polynomial evaluation as an inner product.
3. HPR's selector: after the diagonal-one flip and only for odd primes, every
   `g_1,g_2` instance of the four-term polynomial in Lemma 3.7 is nonzero with
   coefficient in `{+/-1,+/-2}`; the signed level sums of Lemma 3.9 cancel;
   and every non-target family entry supplies a zero prime residue as required
   by Lemma 3.10.

The `Z_6` DGY/Grolmusz-style calibration has `N=6`, theorem-padded dimension
15, and canonical set `{1,3,4}`.  It verifies properties 1 and 2 twice where
applicable.  Property 3 is deliberately marked not applicable because HPR
requires odd primes.  A companion `Z_15` family verifies property 3.

The polynomial generator records both realized degree and the theoretical
`max_i p_i^e_i-1` bound after Boolean multilinearization.  A mismatch is fatal.
This is the degree-accounting guard motivated by withdrawn TR26-044.

## 5. Secondary target after the barrier

A concrete accounting improvement is now implemented in
[`packed_catalyst_lemma.md`](packed_catalyst_lemma.md): jointly encode all
`D=O(d)` ring coordinates as one base-`m` integer.  The arbitrary initial
`ceil(log_2(m^D))`-bit string differs from a valid packed vector by at most one
multiple of `m^D`, so a single persistent quotient bit suffices.  This reduces
HPR's binary catalytic-space term from `O(d log(dm))` to `O(d log m)`.  It is a
real fixed-`m` `Theta(log d)` improvement, but BDL still forces
`d=Omega(ell log ell)`, so it does not reach L.

The most faithful HPR pivot is their own open problem: materialize or manipulate
matching vectors “on the fly” from a representation much shorter than the
explicit `d` coordinates.  A useful subproblem must support, in `O(ell)` live
space, the exact operations HPR uses:

- add a scalar multiple of `u_a` or `v_a` to a logical register;
- compute inner products of a logical register with `v_r`;
- reverse every update exactly; and
- handle the arbitrary initial catalyst, not merely a zero register.

Logspace coordinate generation alone is insufficient: it streams a vector but
does not compress the evolving arbitrary `d`-coordinate register.

A second, narrower experiment is to port the July 2026 optimal-sparsity
decoding polynomials into HPR's catalytic-information-retrieval interface.  For
the same number of prime factors they suggest reducing an exponential server
index to roughly linear size, which could improve HPR's growing-`t` recursion
tradeoff.  This is only a candidate until perfect correctness, deterministic
additive query form, concatenation, reversible answer/reconstruction, state
size, runtime, and polynomial degree are all proved in the white-box HPR model.
It cannot yield L by itself because it leaves `d` unchanged.

## 6. Withdrawn result exclusion

[ECCC TR26-044](https://eccc.weizmann.ac.il/report/2026/044/) was withdrawn on
2026-04-07 because the polynomial degree per subtree was miscalculated, so the
claimed polynomial runtime did not follow.  No theorem or construction in this
project depends on it.

## 7. Growing-modulus caveat

The fixed-`m` barrier must not be misquoted as a barrier to every explicit MV
route. With fixed `t` but growing primes, BDL's general restricted-family bound
only forces the information-theoretic budget `d log m=Omega(ell)`. Such a
family is inside known general upper bounds. It would still not imply L from
HPR Theorem 3.1 as written, because the free-space term `h log m` grows.

The exact construction and recursion-state obstacles, including the rational-
rank barrier when all primes are too large, are recorded in
[`succinct_state_audit.md`](succinct_state_audit.md). Thus the retired search
target is specifically **fixed modulus, `d=O(ell)`**. A growing-modulus search
is legitimate only if it simultaneously accounts for the per-level
reconstruction-state stack.
