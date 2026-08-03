# Research program after the first formalization pass

Date: 1 August 2026.

## Problem framing

The fixed-modulus explicit-family route to TreeEval in L is closed by the
matching-vector upper bound. Dense catalyst packing improves HPR's exact
binary-catalyst parameter but does not change this conclusion. The remaining
program has two coupled unknowns:

1. a uniform canonical MV family with `d log m=O(ell)` in a growing-modulus
   regime; and
2. a same-type recursive reconstruction procedure whose suspended state is
   independent of `log m` and whose branching is independent of `d` and
   `2^ell`.

Either object without the other is insufficient.

## Ranked candidates

| Rank | Candidate | Expected value | Effort | Current evidence |
|---:|---|---|---|---|
| 1 | Nonlocal same-type CRT batch transducer | Highest: directly attacks the `h log m` stack | Very high | Projected components are exact; all local and sequential compositions in the ledger fail |
| 2 | Trace-aware register representation | High: could evade the generic profile barrier | High | The zero-base scan exposes an arbitrary vector-label sequence, so success must exploit recomputation and control, not a small update alphabet |
| 3 | Growing-modulus canonical Gram family | High but only jointly useful with Rank 1 or 2 | Very high | Known upper bounds permit `d log m=Theta(ell)`; the first exceptional finite target `(3,5),d=3,N=10` remains unresolved |
| 4 | Same-type conversion for derivative PIR | Medium: improves server branching if conversion is batched | High | Exact `GF(4)/Z_6` pair kernel works; field-to-exponent conversion currently incurs `d` child replays |
| 5 | Further fixed-modulus MV search | Low | Unbounded | Retired by BDL plus bounded-torsion PFR; retain only calibration instances |

## Experiments and gates

### Milestone 1: two-level symbolic closure harness

Extend the `Z_15` selector verifier from one algebraic level to a symbolic
two-level oracle trace. A candidate passes only if it:

- matches the direct HPR full update for every base inner-product pair and
  canonical delta pattern;
- returns the same `Z_15^d` update type it consumes;
- restores all mask and scratch registers; and
- uses a constant number of child calls independent of the six labels.

**Go:** all identities pass in two independent implementations and the child
call recurrence is `2^O(t h)`. **No-go:** any `N`, `d`, or `log m` factor
multiplies at every level.

### Milestone 2: restricted-trace state/time frontier

Model the zero-base all-active scan and later selector phases as transparent
program prefixes. Measure the exact query-profile rank reachable at each
prefix and the recomputation recurrence when the prefix, rather than the state,
is retained.

**Go:** a checkpoint schedule uses `O(ell+h)` bits and polynomial total replay.
**No-go:** the recurrence contains `d`, `N`, or `N^2` at `Theta(h)` levels.

### Milestone 3: exceptional Gram search — completed no-go

Keep `(p_1,p_2,d,N)=(3,5,3,10)` as the first bounded construction target, but
add symmetry breaking and independently check every witness by explicit
factorization, both MV verifiers, and the Hadamard/rank verifier.

**No-go obtained 2 August 2026:** a complete rank split gives 38 exact
bit-vector cases, all solver-certified UNSAT with reproducible encoding hashes
and no timeout.  The same normalized encoding recovers the `N=9` positive
control.  Therefore the maximum at `(3,5),d=3` is exactly 9; see
`exceptional_gram_n10_no_go.md`.  The zero-one determinant bound lifts this
from `(3,5)` to every pair of distinct odd primes: the dimension-3 maximum is
always exactly 9.

### Milestone 4: derivative-PIR conversion

Search only for conversions that batch an entire field-valued answer vector
into one same-type reversible update. Coordinatewise decoding is excluded by
the existing `d T(h-1)` recurrence.

**Go:** constant child calls per answer vector and exact arbitrary-offset
cancellation. **No-go:** a hidden transcript, initialized accumulator, or
coprime-group homomorphism is assumed.

## Metrics

Every proposed construction is scored on:

1. exact suspended bits per recursive level;
2. child-oracle branching per level;
3. output/query type equality;
4. total distinguishable dynamic profiles;
5. restoration for every initial catalyst;
6. uniform generation space; and
7. a fully substituted runtime recurrence.

The degree and recurrence are expanded before asymptotic simplification. This
is the explicit guard against repeating the error behind withdrawn TR26-044.

## Dependencies and risks

- The packed-catalyst theorem is pending author feedback; none of the primary
  targets depends on treating it as novel.
- The generic profile barrier is not a trace-aware lower bound. Overstating it
  would incorrectly retire the only plausible succinct route.
- The finite-state holonomy result is local-model evidence, not a general CIR
  lower bound.
- A growing-modulus family alone does not improve TreeEval while the recursive
  state remains `h log m`.

## Immediate next implementation

The highest-value next code artifact is the two-level symbolic closure harness.
It provides a single falsifiable interface for nonlocal CRT batching,
projective/scaled update types, and derivative-PIR conversions, and rejects
candidate identities before any asymptotic claim is made.
