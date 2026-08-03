# Joint-packing lemma for HPR's catalytic registers

This note isolates a concrete improvement to the catalytic-space accounting in
HPR Remark 2.3.  It does not improve the matching-vector dimension and does not
put TreeEval in L.

The machine-level proof, including uniform bit evaluation and the logical-swap
invariant, is in `packed_catalyst_formal.md`. Exact consequences for HPR's
corollaries are tabulated in `packed_parameter_consequences.md`.

## Lemma

Suppose a catalytic algorithm treats `D` ring coordinates as a vector in
`Z_m^D`, uses coordinatewise additions by streamed vectors, reads coordinates,
computes inner products modulo `m`, and restores the ring vector at termination.
Assume `D log m` is polynomially bounded in the input length and the added
vectors are generated in space `s`.

Then the ring-vector catalyst can be implemented on an arbitrary-initialized
binary catalytic tape using

`b = ceil(log_2(m^D)) = ceil(D log_2 m)`

catalytic bits and `O(s + log(D log m) + log m)` free bits, with polynomial
arithmetic overhead.  One of those free bits must persist while the catalyst is
live.

Applied jointly to HPR's constant number of `Z_m^d` registers, this changes the
binary-catalyst term from `O(d log(d m))` in Remark 2.3 to `O(d log m)`.
For fixed `m`, that removes a factor `Theta(log d)` from the stated catalyst
encoding overhead.

## Exact state map

Let `M=m^D` and interpret the arbitrary `b`-bit catalyst as an integer
`y in [0,2^b)`.  Minimality of `b` gives `2^b < 2M` unless `M` is itself a
power of two (where `2^b=M`).  Euclidean division therefore gives

`y = q M + x`, where `q in {0,1}` and `x in [0,M)`.

Store `q` in one clean bit and replace the catalyst by `x` (subtract `M` iff
`q=1`).  Interpret

`x = sum_{i=0}^{D-1} x_i m^i`

as the dense base-`m` representation of the ring vector.  After the simulated
algorithm restores `x`, add `qM` and erase `q`.  The exact original `y` is
restored for every initial bit string; there are no invalid encodings.

## Operations and space

- Read coordinate `i` as `floor(x/m^i) mod m`.
- Replace coordinate `x_i` by `x_i+a mod m` by adding
  `((x_i+a mod m)-x_i)m^i` to the packed binary integer in place.
- Stream over the `D` digits to compute an inner product modulo `m`.
- Treat swaps among HPR's constant number of registers as a constant-size
  logical permutation of coordinate blocks.

The required powering, iterated multiplication, and division operations have
DLOGTIME-uniform constant-depth threshold circuits.  A bit of their output can
therefore be evaluated in logarithmic space by recursively evaluating the
constant-depth circuit and counting the inputs to a threshold gate.  See
Hesse–Allender–Barrington,
[*Uniform Constant-Depth Threshold Circuits for Division and Iterated
Multiplication*](https://www.cs.rutgers.edu/~allender/papers/division.pdf).

For a destructive update of digit `i`, first compute and retain the old digit
using `O(log m)` bits.  The signed long addend is then fixed:

`c=((x_i+a mod m)-x_i)m^i`.

Generate the bits of `|c|` on demand using the uniform powering/multiplication
circuit, and make one low-to-high pass over the catalytic bits, overwriting each
bit using one carry or borrow bit.  No bit of the pre-update packed integer is
needed after the head has passed it.  The numerical result remains in
`[0,m^D)` because this operation changes exactly one base-`m` digit.

For a read or inner product, the catalytic tape remains unchanged while the
division circuit is evaluated.  The coordinate index uses `O(log D)` bits, the
modular accumulator and retained digit use `O(log m)` bits, and circuit-gate
addresses/counters use `O(log(D log m))` bits.  Runtime is polynomial in
`D log m`.

Intermediate tape operations need not be logically reversible: catalytic
computation requires restoration of the catalyst at the end, not reversibility
of every transition.  The retained quotient bit makes the normalization map
injective, and the simulated ring algorithm's restoration contract supplies
the final `x` needed to invert it.

## HPR corollary

HPR Theorem 3.1 uses a constant number of `Z_m^d` registers.  Concatenate their
coordinates into one vector with `D=Theta(d)` and apply the lemma once, rather
than validating each coordinate separately as in HPR Remark 2.3.  Swaps of the
three named registers are maintained as a permutation of three logical blocks.
Uniform matching-vector coordinates are already generated within HPR's free
space bound.

Consequently the catalytic-space hypothesis in HPR Theorem 3.1 can be replaced
by `O(d log m)` bits.  Its free space remains
`O(ell+h log m+log d)` and its runtime remains
`poly(2^(ell+h t))`: the added arithmetic is polynomial in `d log m`, which HPR
assumes is polynomially bounded in the TreeEval input size.

For HPR Corollary 3.2 this removes a polynomial factor in the MV dimension from
the exact catalyst bound, although the coarser notation `exp(O(ell^epsilon))`
is unchanged after reparameterizing constants.

## Mechanical check

`src/treeeval_mv/packed_catalyst.py` implements the integer state map,
base-`m` coordinate updates, vector addition, and inner products.
`src/treeeval_mv/bit_tape_packing.py` separately implements normalization and
updates as destructive bit passes with generated long constants.  The test
suite exhausts all 128 possible 7-bit catalysts for `m=3,D=4`, every coordinate,
and every modular delta in both implementations (1,536 transitions apiece), and
also tests the exact-power-of-two case.  This is a finite check of the state
map and destructive transition ordering, not a formalization of the uniform
arithmetic-circuit theorem.

## Consequence and limit

With HPR's Grolmusz/DGY family at fixed `m`, the improved catalyst bound is
`O(d)` rather than `O(d log d)`.  This is a genuine catalytic-space improvement.
It does not evade the BDL+PFR bound `d=Omega(ell log ell)`, so the ordinary
explicit-vector route still cannot reach `O(log n)` total space.

Buhrman–Cleve–Koucký–Loff–Speelman explicitly remark after their Lemma 15 that
stronger high-order-bit compression should reduce their non-power-of-two
simulation to linear catalyst size, but they do not give that compression.
The joint dense encoding above is best viewed as a concrete realization of
that stated target for HPR's product ring. HPR use a coordinatewise offset and
Cook–Pyne use enlarged individual registers with a common shift. Author
feedback is still required on whether this exact realization is already
standard and on its compatibility with HPR's machine accounting; no standalone
novelty claim is made.
