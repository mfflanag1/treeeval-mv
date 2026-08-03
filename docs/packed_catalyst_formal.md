# Formal packed-catalyst simulation

Status: theorem-level draft, audited against the catalytic Turing-machine model
and HPR Theorem 3.1 on 1 August 2026. Author feedback is pending.

## 1. Machine-level statement

Let `A` be a catalytic Turing-machine procedure whose logical catalytic state
is a vector `x=(x_0,...,x_(D-1)) in Z_m^D`. Assume:

1. `A` restores `x` for every initial vector;
2. its access to `x` consists of coordinate reads, coordinate additions,
   streamed inner products modulo `m`, and permutations among a constant number
   of fixed-length coordinate blocks;
3. the coordinate index, streamed update value, and any public vector entry
   can be generated in `s` work bits; and
4. `B=D ceil(log_2 m)` is polynomially bounded in the ordinary input length.

Then `A` can be simulated on an arbitrary-initialized binary catalytic tape
using

`b=ceil(log_2(m^D)) <= B`

catalytic bits and

`O(s+log B+log m)`

work bits, including one bit that persists for the lifetime of the logical
vector. The simulation restores every initial binary catalyst exactly. Its
running time is the running time of `A` times a polynomial in `B`.

The theorem is deliberately restricted to the operations HPR use. It is not a
claim that an arbitrary RAM algorithm receives constant-time random access to
a densely packed non-power-of-two alphabet.

## 2. Total state map

Put `M=m^D` and interpret the `b` catalyst cells, least significant bit first,
as `y in [0,2^b)`. Minimality of `b` gives `2^b<2M`, except when `M` is a power
of two, where `2^b=M`. Hence there is a unique decomposition

`y=qM+x`, where `q in {0,1}` and `0<=x<M`.

Compare `y` with `M`, retain `q` in one work bit, and subtract `M` in place iff
`q=1`. The catalyst now holds `x`, which has the unique dense expansion

`x=sum_(i=0)^(D-1) x_i m^i`, with `x_i in {0,...,m-1}`.

After the logical algorithm restores this vector, add `M` iff `q=1`. The tape
again contains the original `y`. The quotient bit can then be cleared. The map
is total on all `2^b` binary tapes; it does not assume a valid initial ring
encoding and has no unused codewords.

## 3. Coordinate access

The `i`th logical coordinate is

`x_i=floor(x/m^i) mod m`.

Integer division and iterated multiplication have DLOGTIME-uniform `TC^0`
circuits. A bit of a logspace-uniform constant-depth threshold circuit can be
evaluated in logarithmic space: recursively evaluate its constant-depth gate
tree, enumerate the polynomially many inputs of a threshold gate, and retain
only a gate address and a logarithmic counter. Thus a requested bit of `m^i`,
`floor(x/m^i)`, or a product can be generated in `O(log B)` space while the
catalytic tape is scanned without modification.

This is the only use of the Hesse–Allender–Barrington arithmetic theorem. It
does not assume that the whole quotient or power is stored on the work tape.

## 4. Destructive coordinate update

To apply `x_i <- x_i+a mod m`, first compute and retain the old digit `x_i` in
`O(log m)` bits. Let

`x_i'=(x_i+a) mod m` and `c=(x_i'-x_i)m^i`.

The desired packed state is exactly `x+c`. Generate the bits of `|c|` on
demand and make one least-significant-to-most-significant pass over the
catalytic tape, using a single carry bit if `c>=0` and a single borrow bit if
`c<0`. Once a tape bit is overwritten, no later step needs its old value.

There is no underflow or overflow: replacing one valid base-`m` digit by
another leaves the numerical result in `[0,M)`, which is contained in the
`b`-bit tape. Applying the additive inverse repeats the same construction and
restores the previous packed integer.

A streamed inner product retains one `O(log m)` accumulator and obtains each
digit in turn. A permutation among HPR's three logical registers is represented
by one of six block permutations; recursive swaps mutate this single global
permutation and undo it on return. No stack of permutations is needed.

## 5. HPR corollary

HPR's proof uses three registers `x,y,z in Z_m^d`. Concatenate them into one
logical vector with `D=3d` and apply the simulation once. Uniform matching-
vector entries are already generated within HPR's work-space allowance.

The binary catalyst can therefore have length

`ceil(log_2(m^(3d))) = ceil(3d log_2 m) = O(d log m)`.

The persistent quotient bit is global, not one bit per recursion level. Digit
access, inner products, and HPR's final read of the padded output coordinates
fit within

`O(ell+h log m+log(d log m))=O(ell+h log m)`

free space under their polynomial-size hypothesis. Each simulated ring
operation incurs only polynomial overhead in `d log m`, so the theorem's
`poly(2^(ell+h t))` running-time form is unchanged.

Thus the catalyst hypothesis in HPR Theorem 3.1 can be sharpened from
`O(d log(dm))` to `O(d log m)`, subject to confirmation that their intended
machine accounting permits polynomial-time sequential access to the packed
registers.

## 6. Independent finite checks

Two implementations use different state transitions:

- `packed_catalyst.py` performs exact integer normalization and digit updates;
- `bit_tape_packing.py` performs destructive ripple-carry/borrow passes on an
  explicit list of bits.

For `m=3,D=4`, each implementation exhausts all 128 initial seven-bit tapes,
all four coordinates, and all three modular updates: 1,536 forward/inverse
transitions. The exact-power-of-two case is tested separately. These checks
validate the finite state map and overwrite order; the uniform arithmetic
simulation is proved above rather than inferred from testing.

## 7. Scope

- The result improves arbitrary-bit catalytic encoding. A classical clean
  simulation could initialize valid `Z_m` coordinates directly and already
  use `Theta(d log m)` clean bits.
- The result does not improve MV dimension and does not put TreeEval in L.
- BCKLS explicitly anticipated a linear-size high-order-bit compression after
  their Lemma 15. This proof is best treated as a concrete realization for
  HPR's `Z_m` product registers until the authors clarify prior art.

Primary references: HPR Remark 2.3 and Theorem 3.1; BCKLS Lemma 15 and the
paragraph immediately following it; Hesse–Allender–Barrington on uniform
division and iterated multiplication.
