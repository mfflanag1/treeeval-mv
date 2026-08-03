# Stateless four-call closure: exact barrier

Date: 2 August 2026.

## Candidate and why it would have closed

Consider replacing HPR's primewise coefficient filters by a fixed lookup
`phi: Z_m -> Z_m`.  For base inner products `g1,g2` and MV deltas `d1,d2`, its
four shifted child calls have coefficient

`phi(g1*g2) - phi((g1+d1)*g2) - phi(g1*(g2+d2)) + phi((g1+d1)*(g2+d2))`.

If this were one exactly when `d1=d2=1` and zero for every other canonical
CRT-binary delta pair, it would return a full same-type `Z_m^d` update, retain
no selector descriptor across a child call, and compose in 16 leaf calls over
two levels.  Thus it meets the type, state, and recurrence gates of the
two-level closure specification.

## Three-equation obstruction

Let `m` have at least two coprime CRT factors.  Choose a nontrivial CRT
idempotent `e`, and put `f=1-e`; both are non-target canonical deltas.  Set
`g1=0` and `d1=1`.  The mixed difference reduces to
`phi(g2+d2)-phi(g2)`.  Exactness on the following three instances requires

1. `(g2,d2)=(1,e)`: `phi(1+e)-phi(1)=0`;
2. `(g2,d2)=(1,1)`: `phi(2)-phi(1)=1`; and
3. `(g2,d2)=(1+e,f)`: `phi(2)-phi(1+e)=0`.

The first and third equations imply `phi(2)-phi(1)=0`, contradicting the
second over every nontrivial output ring.  For `m=15`, one may take
`{e,f}={10,6}`, producing the explicit cycle `1 -> 11 -> 2` versus `1 -> 2`.

This rules out every stateless scalar reweighting of the existing four HPR
shifts, not merely polynomial or exponent weights.  It does not rule out a
nonlocal transducer with auxiliary state, extra shift patterns, or a changed
recursive type.

## Exact harness

`src/treeeval_mv/fixed_weight_closure.py` implements the mixed-difference
identity and symbolic two-level composition.  The prime-field calibration
`phi(x)=x` over `F_3` passes all 36 one-level contexts and all 1,296 ordered
two-level context pairs with exactly 16 symbolic leaf calls.  The same module
constructs the general CRT contradiction above.
