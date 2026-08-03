# Prior-art log for packed catalytic registers

Checked 2026-08-01.

## Direct antecedent

Buhrman, Cleve, Koucky, Loff, and Speelman, *Computing with a Full Memory:
Catalytic Space* (STOC 2014), Lemma 15, represent non-power-of-two ring
registers using a quadratic number of blocks. Immediately after the proof they
state that stronger compression of the high-order bits could reduce the
catalyst from `(r ceil(log_2 |R|))^2` to `O(r ceil(log_2 |R|))` bits. They do
not provide that compression.

The joint rank encoding implemented here realizes that stated target for a
constant number of registers over `R=Z_m^D`: encode the entire register tuple
as one integer in `[0,m^D)`, normalize an arbitrary `ceil(log_2 m^D)`-bit tape
with one quotient bit, and access base-`m` digits in logspace.

## Later treatment checked

- Henzinger-Pyne-Ragavan (ECCC TR26-022), Remark 2.3, uses coordinatewise
  validity offsets and states `O(d log(dm))` binary catalyst space.
- Cook-Pyne, *Efficient Catalytic Graph Algorithms* (arXiv:2509.06209), uses
  enlarged individual registers plus a common random shift to support fast
  register access. It does not jointly rank-encode the entire tuple.
- Exact web searches for the distinctive BCKLS phrase "stronger compression
  of the high order bits" found copies and a thesis reproduction of the 2014
  remark, but no later construction implementing it.

This is evidence for an omitted implementation, not a novelty proof. The note
to HPR explicitly asks for prior-art and machine-model correction.
