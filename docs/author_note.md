# Questions on the matching-vector TreeEval construction

**Notes for the authors — 1 August 2026**

I wrote exact checkers for the matching-vector and selector identities while
working through the paper. The main question that came out of this is about
the binary representation of the catalytic registers. I have also included
two smaller parameter observations for context.

## 1. Can the catalytic registers be packed jointly?

HPR Remark 2.3 validates `Z_m` coordinates separately and gives
`O(d log(dm))` binary catalytic space. It seems possible to avoid the extra
`log d` factor by packing all `D=Theta(d)` coordinates into one integer.

Let `M=m^D` and `b=ceil(log_2 M)`. Interpret the arbitrary initial `b` catalyst
bits as an integer `y in [0,2^b)`. Since `2^b<2M`, there is a unique

`y=qM+x`, with `q in {0,1}` and `0<=x<M`.

Save `q` in one work bit, subtract `qM` in place, and read
`x=sum_i x_i m^i` as the vector in `Z_m^D`. Once the simulated ring algorithm
has restored `x`, adding `qM` restores the original binary tape exactly.

Coordinate `i` is `floor(x/m^i) mod m`. If it is changed by `a`, retain the
old `O(log m)`-bit digit and add the signed integer

`((x_i+a mod m)-x_i)m^i`

to the packed tape in a low-to-high carry or borrow pass. Bits of the powers,
products, and quotients can be generated in logspace using the uniform TC0
arithmetic of Hesse–Allender–Barrington. Reads and inner products just stream
over the packed digits. HPR use only a constant number of registers, so their
coordinates can be concatenated and register swaps handled by a constant-size
logical permutation.

This appears to use `ceil(D log_2 m)` catalyst bits and
`O(log m+log(D log m))` additional work bits, with polynomial overhead. In
Theorem 3.1 it would replace `O(d log(dm))` by `O(d log m)` without changing
the stated free-space or polynomial-time bounds.

BCKLS remark after Lemma 15 that stronger compression of the high-order bits
should reduce their non-power-of-two simulation to linear catalyst size. I do
not view the target as new; the question is whether the state map above is a
valid realization for HPR's product-ring registers, or whether packed
coordinate access violates a constraint of the machine model.

I checked the state map two ways: an integer implementation and a separate
destructive bit-tape implementation. For `m=3,D=4`, both exhaust all 128
initial seven-bit tapes, every coordinate, and every modular update (1,536
update/inverse transitions in each implementation). The power-of-two boundary
is checked separately.

## 2. The fixed-modulus dimension boundary

Let `N=2^ell`. For constant `m,t`, using Theorem 3.1 as an ordinary clean-space
algorithm would require `d=O(ell)` to reach `O(log n)` total space.
Bhowmick–Dvir–Lovett Theorem 2, together with the now-proved bounded-torsion
PFR theorem of Gowers–Green–Manners–Tao, gives for fixed `m`

`MV(m,d) <= exp(c_m d/log d)`.

Thus `N=2^ell` forces `d=Omega_m(ell log ell)`. This is consistent with the
`d=O(ell log ell)` target in HPR Remark 3.4, which would give polynomial-time
`O(log n log log n)` space, but it rules out reaching L through a fixed-modulus
explicit family alone. My reading is that the remaining L route must use the
vectors through a succinct representation rather than materializing all `d`
coordinates. Is that also how you view the boundary?

## 3. A characteristic-local selector variant

There is a small simplification if the update is projected to one CRT
component `p_k`. At that component, use the shifted product itself instead of
selecting one of its monomials. Its signed mixed difference is

`xy-(x+1)y-x(y+1)+(x+1)(y+1)=1`.

For `i!=k`, keep HPR's `(alpha_i,beta_i)` filters, normalize by their
coefficients, and multiply by the CRT idempotent for `p_k`. Lemma 3.8 can be
run modulo `m/p_k`, so the baseline inner products stored across its four
oracle calls also omit the `p_k` component. The resulting one-level selector
is one at the target and zero elsewhere while suspending
`O(t+log(m/p_k))` bits rather than `O(log m)`.

This identity is confined to its output characteristic. A fixed weight
`phi:F_p->F_q` cannot provide the same mixed difference for distinct primes:
summing over `x in F_p` makes the left side telescope to zero while the right
side is `p!=0` in `F_q`. I checked the positive identity over the `Z_15`
companion family in 7,200 direct and 7,200 independently factorized cases.

I do not see a way to compose the projected updates recursively without
recovering the complementary state, so this is only a one-level observation.
Could it still be useful for an asymmetric modulus, or is this limitation
already captured by the CIR formulation?

## References

1. A. Henzinger, E. Pyne, S. Ragavan, *Catalytic Tree Evaluation From Matching
   Vectors*, ECCC TR26-022, 2026.
2. A. Bhowmick, Z. Dvir, S. Lovett, *New Bounds for Matching Vector Families*,
   SIAM J. Comput. 43(5), 2014.
3. W. T. Gowers, B. Green, F. Manners, T. Tao, *Marton's Conjecture in Abelian
   Groups with Bounded Torsion*, 2024/2025.
4. H. Buhrman, R. Cleve, M. Koucky, B. Loff, F. Speelman, *Computing with a
   Full Memory: Catalytic Space*, STOC 2014.
5. W. Hesse, E. Allender, D. A. M. Barrington, *Uniform Constant-Depth
   Threshold Circuits for Division and Iterated Multiplication*, JCSS 65(4),
   2002.
