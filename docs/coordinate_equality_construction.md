# Coordinate-equality construction baseline

This is a positive, exact, logspace-uniform HPR-convention construction.  It
does not by itself put TreeEval in L, but it identifies a concrete construction
whose algebra is substantially simpler than the general Grolmusz family.

Fix an alphabet size `a >= 2` and `t` distinct odd primes
`p_1,...,p_t`.  Index labels by any `N <= a^t` distinct words
`x=(x_1,...,x_t)` in `[a]^t`.  In CRT component `r`, set

`G_r[x,y] = 1[x_r=y_r]`.

Factor `G_r` using the `a` one-hot vectors.  CRT-combine corresponding
coordinates across the primes.  The resulting vectors lie in `Z_m^a`, where
`m=prod_r p_r`, and satisfy

- diagonal Gram pattern `(1,...,1)`;
- every off-diagonal pattern belongs to `{0,1}^t \ {(1,...,1)}`;
- size `N <= a^t` and dimension `d=a`.

The generator `coordinate_equality_hpr_family` implements the construction.
The test suite checks an instance with `N=27`, `d=3`, and `m=105` using both
independent exact MV verifiers.  A smaller instance is also passed through the
full HPR selector-identity verifier.

## Parameter ledger

Taking the first `t` odd primes gives `log m = Theta(t log t)`.  To provide
`N=2^ell` labels, choose `a=ceil(2^(ell/t))`.  The packed catalyst budget is

`d log m = Theta(2^(ell/t) t log t)`.

At `t=Theta(ell)` this becomes `Theta(ell log ell)`: only a `log ell` factor
above the information target.  However, HPR's current reconstruction has time
`poly(2^(ell+h t))`, so this large-`t` choice is not polynomial time for
general height `h`.  This construction therefore makes the joint target
precise:

1. replace equality kernels by canonical component kernels carrying more than
   `a` labels at rank `a`, while allowing the primes to grow; or
2. exploit the coordinate-equality structure to reconstruct with dependence
   polynomial in `t`, rather than enumerating all `2^t` CRT patterns at every
   recursive level.

The first direction is a finite-rank construction problem.  The second is a
changed-recursion problem, not merely a different encoding of the same HPR
state.

## Uniformity

For a label index, its word is obtained by repeated division by `a`.  A vector
coordinate is the CRT combination of `t` equality bits.  Both operations use
`O(log N + log m)` work space, and selecting the first `2^ell` words preserves
the MV property.  Thus arbitrary powers-of-two family sizes are obtained
without padding or a nonuniform choice.
