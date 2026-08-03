# Author-response playbook

No further message should be sent without explicit permission.

## A. Sound and not already standard

1. Thank the authors and ask how they would like priority and attribution
   handled before public circulation.
2. Offer the formal proof, exact finite checks, and a minimal patch to the HPR
   theorem statement.
3. Prepare a short note centered only on the packing theorem and its precise
   corollaries. Do not market it as TreeEval in L.
4. Ask whether the authors prefer an acknowledgement, a joint revision, or an
   independent note.

## B. Sound but known or intentionally omitted

1. Request the best citation or folklore attribution.
2. Change the project status from “candidate improvement” to “explicit
   derivation of a standard compression.”
3. Keep the implementation as a reproducibility calibration and remove any
   implication of novelty.
4. Continue with the succinct-state and growing-modulus targets.

## C. A machine-model flaw is identified

Classify the objection before changing the theorem:

| Possible objection | Required repair or verdict |
|---|---|
| Packed digit cannot be read in the claimed work space | Give the uniform-circuit bit-evaluation simulation gate by gate; abandon the claim if their model requires constant-time word access |
| Destructive addition needs overwritten high bits | Use the explicit low-to-high carry/borrow invariant and bit-tape test |
| Normalization is not injective | Exhibit `y=qM+x`, retain `q`, and show exact restoration |
| The quotient bit is charged once per recursion level | Show normalization occurs once around the three global HPR registers; `q` is one global bit |
| Register swaps require physical movement | Use the six-state logical block permutation and prove nested swaps restore it |
| Arithmetic overhead breaks polynomial time | Substitute the `d log m <= poly(2^(h+ell))` hypothesis into every simulated operation |
| A hidden extra catalytic register is needed | Reject the implementation unless it fits the same packed `D=3d` state |

Add a regression test reproducing any finite counterexample supplied by the
authors before attempting a repair.

## D. Only a native-ring version is accepted

Separate two statements:

- native `Z_m` catalytic cells need `Theta(d log m)` information bits;
- the proposed theorem handles an arbitrary binary tape with the same
  asymptotic size.

If the binary lift is rejected but the native-ring accounting is accepted,
withdraw the claimed improvement to Remark 2.3 and retain only the native-ring
observation.

## E. No response

After a reasonable interval, prepare one concise follow-up containing only the
formal theorem and the single machine-model question. Do not send it without
new authorization. In parallel, seek an independent complexity-theory review
of the proof rather than treating silence as confirmation.
