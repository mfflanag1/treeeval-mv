# Harvest closeout — 21 August 2026

## Decision

Park the project as an open-ended TreeEval-in-L campaign. Preserve it as a
reproducible matching-vector / catalytic-computing audit and reopen it only on
a concrete new mechanism or external expert feedback.

The project has strong verification infrastructure and several bounded results,
but no scaling construction or recursively polynomial-time vector type. The
two remaining requirements are independent:

1. a growing-modulus canonical matching-vector family with
   `d log m = O(ell)`; and
2. a recursively closed reconstruction with constant effective branching and
   no `h log m` suspended state.

The coordinate-equality family is an exact baseline, not a solution: it loses
a `log ell` factor in packed space and its direct HPR reconstruction has
superpolynomial recursion when the number of prime factors grows. The only
new finite-geometry candidate tested so far is strictly worse than that
baseline.

## Frontier refresh

Henzinger--Pyne--Ragavan's result remains the current positive TreeEval result
used here: logarithmic free space, polynomial time, and subpolynomial
catalytic space from uniform matching vectors. Its stated open direction is
better matching-vector constructions / a more succinct materialization route.

On 10 August, Aggarwal and Obremski posted ECCC TR26-141, *Subexponential Upper
Bounds for 3-Restricted Matching Vector Families*. It proves, for every fixed
modulus, `MV(m,n,3) <= exp(O_m(sqrt(n log n)))`. This is relevant evidence
that the matching-vector frontier is active and that fixed-restriction routes
are tightening. It does **not** directly rule out this project's HPR target:
the canonical CRT family with `t` prime factors has up to `2^t` inner-product
patterns, and the desired escape regime has growing modulus.

The Asadi--Cleve catalytic-pebbling manuscript remains withdrawn. ECCC states
that its polynomial-degree calculation for each subtree was erroneous and its
claimed polynomial runtime therefore does not follow. This directly validates
the project's practice of expanding recurrences before accepting an
asymptotic claim.

## Author-contact gate

The HPR packing note was sent on 1 August. A mailbox search on 21 August found
the sent thread and no reply from the three recipients. Therefore the
`O(d log m)` packing observation remains unaudited externally and must not be
promoted as a novel theorem.

## What is retained

- two independent exact MV verifiers;
- an exact Gram/rank characterization and bounded UNSAT ledger;
- a uniform coordinate-equality growing-modulus baseline;
- an exact implicit-vector and reversal calibration;
- a packed-catalyst implementation and formal accounting; and
- explicit no-go and recurrence gates that distinguish valid finite algebra
  from a TreeEval consequence.

## Reopen conditions

Resume only if at least one of the following occurs:

1. the HPR authors validate the packing lemma as sound and substantively new;
2. a candidate kernel provably beats equality at multiple growing scales with
   `log N = Omega(d log m)`;
3. a two-level symbolic closure achieves constant effective child branching;
   or
4. a collaborator supplies a specific algebraic or recursive ansatz not
   already excluded by the existing ledger.

Absent such a trigger, more search would be infrastructure growth rather than
a proportionate route to the breakthrough.

## Sources checked

- Henzinger, Pyne, Ragavan, ECCC TR26-022 / arXiv:2602.14320.
- Aggarwal, Obremski, ECCC TR26-141, revision 2, 11 August 2026.
- Asadi, Cleve, ECCC TR26-044 withdrawal notice, 7 April 2026.
