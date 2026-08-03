# Exact parameter consequences of joint packing

Status: checked 1 August 2026. These are representation improvements, not a
TreeEval-in-L claim.

Write

`A_t(ell)=ell^(1/t) (log ell)^(1-1/t)`.

| Regime | HPR binary catalyst accounting | Joint packing | Consequence |
|---|---:|---:|---|
| Theorem 3.1, arbitrary uniform MV family | `O(d(log d+log m))` | `O(d log m)` | Free space and `poly(2^(ell+h t))` time unchanged |
| Corollary 3.2, fixed `t,m`, `d=exp(O(A_t(ell)))` | `O(d log d)` | `O(d)` | Removes the representation-level `log d` factor; the published coarse `exp(O(ell^epsilon))` class is unchanged |
| Corollary 3.3, first `t` odd primes, `t=sqrt(log ell-(log log ell)/2+O(1))` | `O(d(log d+log m))` | `O(d sqrt(log ell) log log ell)` | Uses `log m=Theta(t log t)`; both bounds remain `exp(exp(O(sqrt(log ell))))` in coarse notation |
| Hypothetical fixed-`m` family, `d=O(ell log ell)` | `O(ell (log ell)^2)` as an arbitrary-bit catalyst | `O(ell log ell)` | Matches the information stored in the explicit ring vector |
| Information-tight growing modulus, `d log m=O(ell)` | may retain an extra `d log d` term | `O(ell)` | Catalyst size would be compatible with L, but HPR's suspended `h log m` free state and the missing MV family remain |

## Clean-space versus catalytic-space accounting

HPR Remark 3.4 already obtains polynomial-time `O(log n log log n)` classical
space from a hypothetical fixed-modulus `d=O(ell log ell)` family. There is no
contradiction with the fourth row: a classical simulation chooses valid zero
ring registers and pays `Theta(d log m)` bits directly. The larger
`O(d log(dm))` term is needed only when every initial binary catalyst string
must be accepted and restored.

## What changes in the published corollaries

The standard big-Oh headlines of Corollaries 3.2 and 3.3 do not change because
their MV dimensions dominate logarithmic representation factors. A sharpened
version should therefore state the exact catalyst dependence

`C_HPR(d,m)=O(d log m)`

before substituting the construction. The gain is material for fine-grained
catalyst accounting and for any future growing-modulus family satisfying
`d log m=Theta(ell)`.

## What does not change

1. The fixed-modulus BDL plus bounded-torsion PFR barrier still forces
   `d=Omega_m(ell log ell)`.
2. Packing supplies no succinct dynamic representation of a register smaller
   than its information content.
3. When `m` grows, HPR's current recursive selector still suspends
   `Theta(log m)` reconstruction data per level.
