# Succinct operational-state audit

Status: 2026-08-01. This note records the shortcuts that were tested and the
exact point at which each one fails. None is a TreeEval-in-L claim.

## 1. The short descriptor is not a fast representation

An HPR execution prefix has a short control description. At every ancestor
oracle call, the caller retains only a server/phase symbol of `O(t)` bits and
the node path takes `O(h)` bits. At the deepest active level, the public
truth-table counters `r,s` use `O(ell)` bits. Thus the control prefix is already
within `O(ell+h t)` space.

That prefix does not provide random access to the three vector registers. A
single one-level call directly performs as many as

`2^(2t) * 2^(2ell)`

updates of `z` by public vectors `w_{f(r,s)}`. The cancellations proving that
the final sum is one matching vector happen only after the full selector sum.
A sparse list of uncombined terms is therefore exponentially long even before
recursion.

## 2. Coordinate replay is superpolynomial

One can recompute a requested register coordinate by replaying the short
program prefix and recursively materializing each child contribution only for
that coordinate. This saves space but multiplies work across levels.

Two natural implementations have recurrences of the following form:

- coordinate replay through an explicit inner-product scan:
  `T(h) >= d T(h-1)`;
- stateless Gram replay, which recomputes the child feature separately for
  each public candidate pair: `T(h) >= 2^(2ell) T(h-1)`.

For `h=Theta(ell)`, the first is superpolynomial for every growing `d`, and the
second is `2^(Theta(ell^2))`. Checkpointing interpolates between these two
failures: polynomial replay needs checkpoints at constantly spaced levels,
which stores `Theta(h ell)` label bits; `O(ell+h)` space leaves a block of
`Theta(h)` levels and restores the superpolynomial replay.

This is a no-go for these replay strategies, not a lower bound against every
operational representation.

## 3. Direct label/XOR registers recreate the same tradeoff

If a logical register stores an `ell`-bit label and the child oracle XORs the
unknown child value into it, one bit of the child value is easy to recover by
reading the register before and after a balanced oracle call. Recovering both
child labels and retaining them while recursing costs `Theta(h ell)` bits.

Scanning the truth table and recomputing only the currently needed equality
bit uses `O(ell+h)` space, but it makes `2^(2ell)` recursive calls per level and
has the second recurrence above. Matching vectors are precisely a batched
equality sketch that avoids this recomputation; the explicit sketch is what
costs `d log m` bits.

## 4. Growing modulus: an inside-the-upper-bound target, but two obstacles

The fixed-modulus BDL/PFR barrier does not apply when the prime factors grow
with `ell`. For a canonical family with a fixed number `t` of prime factors,
BDL Theorem 1 gives only

`N <= C_t m^(d/2)`.

Consequently `N=2^ell` forces `d log m >= 2ell-O_t(1)`. This permits, at least
at the level of known upper bounds, an information-tight target
`d log m=Theta(ell)`.

There is a further elementary large-prime barrier. Work directly in HPR's
diagonal-one convention and let `G_r` be the zero-one Gram matrix modulo
`p_r`. It has rank at most `d` over `F_{p_r}`, and

`G_1 circ ... circ G_t = I_N`.

The sharp zero-one determinant inequality says that an `s` by `s` minor has
absolute value at most `(s+1)^((s+1)/2)/2^s`. Therefore, if the exact integer
comparison

`(p_r 2^(d+1))^2 > (d+2)^(d+2)`

holds, every `(d+1)` minor that vanishes modulo `p_r` also vanishes over the
integers, so `rank_Q(G_r)<=d`. Rational Hadamard rank then gives

`N <= d^t`.

Starting from a customary dimension-`d_0` family and applying HPR's affine
flip gives `d=d_0+1`, recovering the earlier `(d_0+1)^t` form. Thus an
information-tight growing-modulus construction must keep at least one prime
below `exp(O(d log d))`. For balanced primes this pushes the first plausible
regime to roughly

`d = Omega(sqrt(ell/log ell)),  log m = O(sqrt(ell log ell))`.

No canonical construction near this frontier was found in the large-modulus
MV literature; the known near-`m^(d/2)` constructions allow a large avoiding
set and do not satisfy HPR's constant CRT cube.

The Gram characterization is also an exact finite search object. It proves
`N_max=4` for `d=2,t=2` and `N_max=8` for `d=2,t=3`, attained by the
coordinate-equality product construction. The first odd-characteristic rank
exception is `d=3,p=3`; the smallest capacity-improving target is a
`(3,5),d=3,N=10` Gram system. Its bounded bit-vector search currently times
out, so it remains unresolved rather than entering the no-go ledger.

Even such a family would not plug directly into HPR. Their recursion stores
`Theta(log m)` reconstruction data per level, giving `h log m` free space.
The tempting fix is to retain only constant-size CRT bit patterns and
recompute `g,alpha,beta` after every child call. It fails in its naive form:
the routine computing `g=<x,v_a>` itself must preserve a `log m`-bit baseline
inner product across recursive oracle calls. Bit-serial recomputation removes
the baseline but introduces `log m`-dependent branching per level, yielding
superpolynomial time when `h=Theta(ell)` and `m` grows.

The growing-modulus escape therefore requires **both** a near-information-
tight canonical family and a genuinely transparent reconstruction-state stack.

## 5. Exact remaining breakthrough object

Before looking for that object, there is an information-theoretic barrier to
generic dynamic compression. If `G_i` is the component Gram matrix modulo
`p_i`, arbitrary MV updates from zero generate exactly
`p_i^rank(G_i)` distinct complete query profiles in that component. CRT makes
the components independent. Hence any exact state representation supporting
arbitrary updates and all Gram queries—even a nonlinear one—needs at least

`ceil(log_2(prod_i p_i^rank(G_i)))`

bits. See `dynamic_profile_barrier.md` for the proof and an exhaustive `Z_15`
calibration. A successful representation must therefore exploit the restricted
HPR trace, pay for recomputation, or change the interface; a short generic
dynamic sketch cannot suffice.

A useful succinct operational representation must support a batched Gram
transducer, not just a short static label:

1. accept an implicitly evaluated child value and a logical register state;
2. enumerate all values `<X,v_r>` needed by one truth-table scan in total
   polynomial time;
3. support balanced `u_a/v_a` updates and exact reversal;
4. use `O(ell+h)` live bits, including all suspended recursive state; and
5. avoid a `d`, `2^ell`, or `log m` multiplicative factor at every tree level.

The derivative-PIR lift in `cir_pir_audit.md` improves branching once such
state exists, but does not supply this transducer.

## 6. Characteristic-local state can be removed, but not globally

The follow-up in `projected_selector_lemma.md` sharpens the reconstruction
obstacle. For an output restricted to one CRT component `p_k`, the selector
can replace coefficient extraction at `p_k` by the mixed finite difference of
`xy`. It therefore suspends only

`O(t + log(m/p_k))`

bits across oracle calls. An exact `Z_15` verifier checks 7,200 global and
7,200 factorized identities.

The same trick cannot carry a source component of characteristic `p` into an
output component of distinct characteristic `q`: summing the required fixed
mixed-difference identity over `x in F_p` gives `0=p` in `F_q`. Thus the
remaining transducer is specifically **cross-characteristic**. It must close
the separately correct projected updates under recursion without reinstating
the complementary-prime state at every level.

The natural stateful local repair is now closed as well. A complete source
cycle accumulates target-field holonomy `p`. If a reversible persistent state
returns after `L` such cycles, exact restoration requires `L p=0 mod q`, hence
`q|L` for distinct primes. The transducer therefore needs at least `q` states,
or `ceil(log_2 q)` bits. This lower bound is exhausted independently for
`F_3 -> F_5` and `F_5 -> F_3` in `holonomy_barrier.py`. A surviving approach
must batch nonlocally across candidates/components or change the recursive
interface; a small local automaton cannot supply the missing closure.

Primary bounds used here: [BDL](https://arxiv.org/abs/1204.1367) and
[Dvir-Hu large-modulus MV bounds](https://arxiv.org/abs/1304.4819).
