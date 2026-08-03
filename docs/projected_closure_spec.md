# Recursive closure specification for projected selectors

Status: exact target and accounting ledger, 1 August 2026. No recursively
closed construction is claimed.

## 1. The verified primitive

For one output prime `p_k`, the projected selector performs the correct update
in that CRT component. Run HPR Lemma 3.8 modulo `m/p_k`; the baseline inner
products stored across its oracle calls then omit the output component. The
primitive suspends

`O(t+log(m/p_k))`

bits and passes the exhaustive `Z_15` identity checks.

Call this primitive `ProjectedUpdate(k)`. Summing `ProjectedUpdate(k)` over all
`k` gives the desired full `Z_m^d` update at one level.

## 2. Why sequential composition does not close recursively

At a recursive oracle call, the parent routine's suspended state remains live
through the whole child computation. A full update runs every output component,
so some recursion path enters the component with the largest complement. Its
per-level state is

`max_k ceil(log_2(m/p_k)) = ceil(log_2(m/p_min))` up to rounding.

Thus the obvious sequential implementation uses

`Omega(h log(m/p_min))`

stack bits. For balanced or moderately unbalanced primes this is still
`Theta(h log m)`. Component ordering does not help: space is a worst-path
measure, and every component is eventually executed.

For `m=15`, the exact complement costs are three bits for output prime 3 and
two bits for output prime 5. At depth seven, the ordinary four-bit modulus
ledger uses 28 stack bits and the sequential projected ledger still uses 21.

## 3. Why candidate-outer composition loses polynomial time

One can avoid retaining all component filters by fixing a candidate pair
`(r,s)`, cycling through the required child shifts, and deciding its complete
coefficient before touching the output register. This makes at least one
batch of child-oracle calls per candidate pair. With `N=2^ell` labels, the
recurrence has branching at least `N^2` per level and therefore at least

`N^(2h)=2^(2 ell h)`

oracle leaves. For `h=Theta(ell)` this is superpolynomial. The exact exponent
ledger is implemented in `projected_closure_accounting.py`.

## 4. Required object

A recursively closed projected transducer must satisfy all of the following:

1. **Full type:** implement an update in `Z_m^d`, not only one CRT component.
2. **Same recursive interface:** its child queries and returned update have the
   same matching-vector type.
3. **Small suspended state:** retain `O(t)` or `O(1)` bits, independent of
   `log m`, across every child-oracle call.
4. **Batched candidates:** make only `2^O(t)` child-oracle calls per level; no
   `2^ell`, `d`, or `log m` factor may multiply at every depth.
5. **Exactness:** work for every catalyst state and every TreeEval truth table,
   and restore all catalytic storage.
6. **Uniformity:** generate every operation in `O(ell+log m)` space and keep
   the total recurrence polynomial in the input size.

These conditions are a concrete search specification. A proposal that omits
any one of them does not plug into HPR Theorem 3.1.

## 5. No-go ledger around the target

| Attempt | Exact failure |
|---|---|
| Sequential characteristic projections | worst complementary state remains live at every depth |
| Candidate-outer reconstruction | `2^(2 ell h)` oracle-leaf lower bound |
| Fixed cross-characteristic weight | telescoping gives `0=p` in the output field |
| Smaller universal coefficient menu | every menu hitting all base pairs has `Omega(p)` exponents, hence needs `Omega(log p)` selector state |
| Small local reversible repair | pays the target-field holonomy in the stated local model |
| Generic dynamic MV sketch | needs `prod_i p_i^rank(G_i)` distinguishable answer-profile states |
| Derivative-PIR field answer | output type differs from the exponent-vector query and coordinate replay contributes a `d` factor per level |

## 6. Surviving directions

The remaining possibilities are genuinely nonlocal:

- batch Boolean zero-pattern tests across CRT components in a common
  same-type algebra;
- alter the recursive type so projected answers compose without converting
  one coordinate at a time;
- exploit the restricted HPR trace with a polynomial-time replay/checkpoint
  scheme that evades the generic profile barrier; or
- pair an information-tight growing-modulus MV family with a reconstruction
  procedure whose suspended state is transparent across recursion.

The first finite test for any proposal should instantiate `Z_15`, compare its
full update with the existing direct HPR selector on all base inner products
and canonical delta patterns, and then measure child-oracle calls under two
levels of symbolic recursion.
