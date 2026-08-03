# Cross-characteristic selector holonomy

Status: exact lower bound in a natural reversible finite-state transducer
model, verified 2026-08-01. This is a barrier for the obvious stateful repair
of the projected selector, not a lower bound for arbitrary TreeEval
algorithms.

## 1. From telescoping impossibility to a state lower bound

Let `p` and `q` be distinct primes. A target-field weight on the source torus
would have to satisfy

`D_x D_y W(x,y) = 1 in F_q`

for all `x,y in F_p`. Summing over one complete `x` cycle makes the left side
zero for a single-valued potential and the right side `p != 0 in F_q`. This is
the fixed-weight impossibility already used by the projected-selector lemma.

Now allow a reversible finite state to carry the failure to telescope. One
complete source cycle has returned `x` to its starting value but has accumulated
the target offset `p`. If the induced return map on the persistent state has a
cycle of length `L`, exact restoration after that state cycle requires

`L p = 0 in F_q`.

Since `p` is invertible modulo `q`, `q` divides `L`. Thus every reachable
reversible state cycle that works for an arbitrary initial state has at least
`q` states, requiring at least `ceil(log_2 q)` persistent bits. If the state is
required to restore after a single source cycle, no protocol exists at all.

More generally, for source modulus `a`, target modulus `b`, and constant cell
increment `c`, the exact minimum cycle length is

`b / gcd(a c, b)`.

The bound is tight at the macro level: a cycle of that length records the
holonomy and returns after the required number of source circuits. This does
not assert that every such macro cycle implements the full HPR selector; it is
only a necessary state cost.

## 2. Consequence for the projected HPR selector

The characteristic-local selector saves the coefficient-extraction state at
its output prime because its mixed difference is evaluated in the same
characteristic. Trying to replace another prime's Boolean coefficient filter
by a local stateful mixed-difference lift into output field `F_q` necessarily
reintroduces `log q` persistent bits. For a dominant output prime, this erases
the hoped-for state advantage.

This closes the simplest recursively closed transducer model:

1. state-independent cross-characteristic weights are impossible by
   telescoping; and
2. reversible finite-state repairs need a full target-field-sized state cycle.

An escape must therefore use nonlocal batching across candidates or
components, alter the recursion interface, or exploit global register
structure beyond a local cross-characteristic transducer.

## 3. Verification

`src/treeeval_mv/holonomy_barrier.py` computes the exact divisibility bound and
independently exhausts every state permutation up to the claimed minimum for
small prime pairs. For `F_3 -> F_5`, state counts 1 through 4 fail and 5
succeeds; for `F_5 -> F_3`, counts 1 and 2 fail and 3 succeeds. The construction
and exhaustion use integer arithmetic only.
