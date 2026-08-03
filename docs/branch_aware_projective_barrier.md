# Branch-aware and projective selector barrier

Date: 2 August 2026.

## 1. Stronger stateless interface

Allow a selector weight to be an arbitrary function of:

- every CRT branch bit;
- both current inner products separately, not merely their product; and
- the full residue in `Z_m`.

The inner products are computed only during the nonrecursive truth-table scan,
so this interface would retain only `O(t)` branch bits across child calls.  It
strictly contains the fixed four-call lookup ruled out earlier.

## 2. Exact normalization is impossible

The obstruction does not depend on HPR's CRT-cube shifts.  Let `S` be any
finite collection of additive shifts and let `A_s` be an arbitrary function
for shift `s`.  Sum the all-component target equation around the first prime
cycle.  In each `s` term, reindex the cycle by the first CRT coordinate of
`s`.  The resulting left side is the sum around that cycle of the
single-second-component-delta equation, hence zero.  The target right side is
`p_1`, which is nonzero in `Z_m` because another coprime factor divides `m`.
A two-input selector restricts to this one-input estimator after fixing its
second argument to a target edge.

Thus no exact stateless selector exists even with arbitrary branch-dependent
access to both live inner products and an arbitrary finite additive child-query
pattern.  An independent exhaustive `Z_15` sweep checked every one of the
`2^14` shift sets containing zero and found no consistent system.

## 3. Why projective scaling almost works

For any odd cycle and odd modulus there is a cyclic potential whose adjacent
differences are all units: use increments `(1,1,-2)` followed by pairs
`(1,-1)`.  A tensor product of these potentials gives zero on every non-target
CRT pattern and a variable unit on an exact unit target step.  This supplies a
positive one-level projective calibration.

Recursive closure requires the next invocation to accept every scale produced
by the previous one.  Any finite multiplicatively closed set of units
containing one is a subgroup.  `Z_15^*` has exactly eight such subgroups.

An exact finite synthesis was run for all eight.  For each candidate scale
group `K`, the homogeneous off-target equations were reduced independently
over `F_3` and `F_5`; target responses were constrained to the corresponding
residues of `K`.  Every group is UNSAT.  The full unit group fails modulo 3;
the smaller candidates that survive nonvanishing fail their required mod-5
target residues.  Hence a scalar projective fixed point does not exist at the
smallest HPR modulus.

This last statement is a bounded computational no-go over `Z_15`, unlike the
general exact-normalization theorem.  It does not rule out a changed vector
type, non-scalar state, or a selector with a different child-query pattern.

## 4. Why one shared arbitrary register is not clean state

Suppose a catalytic register ranges over a finite set `C`.  For every value
`beta`, a reversible attempt to encode `beta` is a permutation `E_beta:C->C`.
Every permutation has image `C`.  Therefore a decoder seeing only
`E_beta(c)`, for arbitrary unknown `c`, cannot recover `beta`: the same encoded
state is reachable for every `beta`.  Another clean value or a correlated
baseline must remain live, recreating the resource the encoding was intended
to remove.

Consequently a single `log m`-bit catalyst threaded through all recursion
levels does not by itself replace HPR's clean `beta` continuation.  A useful
catalytic-state construction must exploit additional algebraic action on the
vector registers, not merely treat arbitrary bits as reusable clean storage.
