# TreeEval matching-vector audit

This directory is a reproducible audit of the matching-vector route in
Henzinger–Pyne–Ragavan (HPR), *Catalytic Tree Evaluation From Matching
Vectors*.  It combines exact construction experiments with a ledger that
retires parameter regimes only after their obstruction is proved.  The
ordinary fixed-modulus explicit-vector route to classical logspace is outside
the now-unconditional matching-vector upper bound, so current construction
work targets growing moduli and a changed recursive vector type.

The construction-first branch now includes an exact
[`coordinate-equality` family](docs/coordinate_equality_construction.md): a
uniform growing-modulus canonical family with `N=a^t` and `d=a`.  Its parameter
ledger isolates a `log t` storage loss and the generic HPR `2^t` reconstruction
loss as the next two constructive targets.

The first end-to-end construction/type/verification pass is recorded in
[`construction_pipeline_2026-08-02.md`](docs/construction_pipeline_2026-08-02.md).
Its structured-kernel factorizer and implicit-vector prototype are executable;
the report keeps polynomial-time composition as a failed gate, not a claimed
consequence.

The campaign is intentionally parked after a dated frontier refresh in
[`harvest_closeout_2026-08-21.md`](docs/harvest_closeout_2026-08-21.md). It
records the new fixed-modulus 3-restricted upper bound, the unanswered author
query, and the concrete gates required before this line should be reopened.

The consolidated outcome and surviving target are in
[`docs/breakthrough_report.md`](docs/breakthrough_report.md).

The code provides:

- a primary exact modular-arithmetic MV verifier;
- a separately implemented primewise/CRT verifier;
- verification of DGY's polynomial-matching representation;
- an exhaustive checker for the selector identity in HPR Lemmas 3.7–3.10;
- an exact characteristic-local projected selector, plus a finite independent
  check of the cross-characteristic obstruction to making it state-free in
  every CRT component;
- an exact binary-Gram characterization, a bit-vector factor search, and an
  independent modular-rank/Hadamard verifier for asymmetric growing moduli;
- a reversible finite-state holonomy lower bound showing that the natural
  stateful cross-characteristic repair needs a full target-field-sized state;
- a materialized DGY/Grolmusz-style calibration over `Z_6`;
- a companion calibration over `Z_15`, because HPR's selector requires odd
  primes and therefore cannot be calibrated over `Z_6`;
- exact structured set-system and unstructured finite searches, retained as
  calibration/no-go data rather than as a route to `d=O(ell)`;
- an independent SMT encoding for each logged unstructured no-go threshold;
- a jointly packed arbitrary-catalyst representation reducing HPR's binary
  encoding overhead from `O(d log(dm))` to `O(d log m)`, with an exhaustive
  finite state-map check;
- an information-theoretic lower bound for generic dynamic MV sketches,
  calibrated by direct profile enumeration; and
- an exact coefficient-menu search showing why a smaller universal set of
  HPR selector monomials still costs logarithmic state.

The mathematical parameter map and search specification are in
[`docs/parameter_barrier.md`](docs/parameter_barrier.md).  Generated results
are in `results/calibration.json`, `results/no_go_ledger.json`, and
`results/gram_search_ledger.json`.
The secondary packing result is stated separately in
[`docs/packed_catalyst_lemma.md`](docs/packed_catalyst_lemma.md).
Its machine-level proof and exact corollary accounting are in
[`docs/packed_catalyst_formal.md`](docs/packed_catalyst_formal.md) and
[`docs/packed_parameter_consequences.md`](docs/packed_parameter_consequences.md).
Possible author-feedback outcomes are handled in
[`docs/author_response_playbook.md`](docs/author_response_playbook.md).
The post-barrier search policy and operational target are in
[`docs/search_retarget.md`](docs/search_retarget.md).
The exact interface audit, finite derivative-kernel calibration, and the
field-to-exponent recursion no-go are in
[`docs/cir_pir_audit.md`](docs/cir_pir_audit.md).
The operational replay, XOR-label, and growing-modulus escape attempts are
audited in [`docs/succinct_state_audit.md`](docs/succinct_state_audit.md).
The sharper projected-selector lemma and its recursive closure gap are in
[`docs/projected_selector_lemma.md`](docs/projected_selector_lemma.md).
The finite-state obstruction to repairing that gap locally is in
[`docs/cross_characteristic_holonomy.md`](docs/cross_characteristic_holonomy.md).
The generic dynamic-state lower bound is in
[`docs/dynamic_profile_barrier.md`](docs/dynamic_profile_barrier.md). The exact
recursive closure specification and coefficient-menu obstruction are in
[`docs/projected_closure_spec.md`](docs/projected_closure_spec.md) and
[`docs/selector_hitting_set_barrier.md`](docs/selector_hitting_set_barrier.md).
The ranked follow-up program and go/no-go gates are in
[`docs/research_program_after_audit.md`](docs/research_program_after_audit.md).

## Reproduce

Python 3.11+ and a `z3` executable are sufficient; the package itself has no
third-party Python dependency.

```sh
make test
make calibrate
make search
PYTHONPATH=src python3 scripts/run_gram_search.py
```

`make search` exhausts the full pair graph for `(m,d)=(6,1)` and `(6,2)`,
then asks Z3 independently whether a clique one larger than the computed
maximum exists.  It also exhausts the 6-uniform restricted-intersection class
on ground sets of size 6 through 9.  These bounded results test the verifiers
and construction templates.  They are not evidence that the asymptotically
impossible fixed-modulus `d=O(ell)` target might still work.

The Gram search works directly with finite-field factors of binary component
Gram matrices. It proves `N_max=4` at `(p_1,p_2,d)=(3,5,2)` both by a rational
rank bound and an independent `N=4` SAT / `N=5` UNSAT pair. Its `d=3,N=10`
timeout is explicitly logged as unresolved, not as a bounded no-go.

## Claim discipline

- “Exhaustive” always names a finite candidate universe in the JSON ledger.
- Every witness is checked by both MV verifiers.
- The `Z_6` object calibrates the customary diagonal-zero MV and polynomial
  matching properties.  It does **not** calibrate HPR's selector identity.
- The July 2026 sparse-decoding-polynomial result changes server/branching
  parameters, not the MV dimension lower bound.  No TreeEval consequence is
  claimed without a separate white-box catalytic-information-retrieval proof.
- Withdrawn ECCC TR26-044 is not used.
- The reviewed author email was sent on 1 August 2026. No author confirmation
  of the packed-catalyst theorem is assumed.
