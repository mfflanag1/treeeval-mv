# Dynamic matching-vector profile barrier

Status: exact theorem and finite calibration, 1 August 2026. This is a barrier
for generic dynamic compression, not for trace-specific TreeEval algorithms.

## Theorem

Let `(U,V)` be an HPR canonical MV family over squarefree
`m=prod_(i=1)^t p_i`, with component Gram matrices

`G_i=U_i V_i^T over F_(p_i)`.

Consider any state representation that starts from the zero vector, supports
arbitrary updates

`X <- X + gamma u_a`, with `gamma in Z_m`,

and from its current state answers every query `<X,v_r> mod m` exactly. The
representation has at least

`prod_i p_i^(rank(G_i))`

distinct states and therefore needs at least

`ceil(log_2(prod_i p_i^(rank(G_i))))`

bits in the worst case. The encoding and query algorithms need not be linear.

## Proof

Fix a prime component `p_i`. From zero, coefficient vector
`c in F_(p_i)^N` reaches `X=c^T U_i`. Its complete vector of Gram-query answers
is

`(<X,v_1>,...,<X,v_N>) = c^T G_i`.

As `c` varies, these profiles form the row space of `G_i`, which has exactly
`p_i^(rank(G_i))` elements. States with different profiles must be distinct,
or some query would receive two different correct answers from the same state.

CRT permits the update coefficient of every `u_a` to be chosen independently
in every prime component. The joint profile set is therefore the Cartesian
product of the component profile sets, and its cardinality is the stated
product.

## Consequence for the succinct HPR target

A context-independent dynamic data structure cannot compress the evolving MV
register below its true component Gram ranks while retaining arbitrary updates
and all exact queries. In particular, replacing an explicit register by a
short Gram oracle is insufficient: the oracle describes a static row, whereas
the dynamic state must distinguish every reachable answer profile.

An L algorithm must evade the theorem in one of three explicit ways:

1. support only the restricted update traces produced by the recursive HPR
   control flow;
2. recompute query answers from that trace instead of storing a generic state;
   or
3. change the recursive interface so that it does not expose arbitrary MV
   updates and all Gram queries.

The theorem does not rule out any of these routes. It does rule out lossless
generic sketches, including nonlinear ones, whose advertised size is below
the answer-profile entropy.

The HPR trace is not narrow merely at the level of its update alphabet. With
both logical masks zero and the all-zero shift, every component selector is
`X-1`, so `beta_i=0` and all `N^2` candidate pairs are active. The arbitrary
truth table can send an arbitrary sequence of matching-vector labels into the
evolving output register during that scan. This does not by itself prove the
full state lower bound for a trace-aware algorithm—the input and scan counter
can be used for recomputation—but it shows why “the trace only uses a few
updates” is not a viable escape.

## Finite calibration

For the coordinate-equality HPR family over `Z_15` with alphabet size two,
`N=4`, `d=2`, and both component Gram ranks are two. Exhausting all coefficient
vectors gives:

- `3^2=9` profiles modulo 3;
- `5^2=25` profiles modulo 5; and
- `225` joint profiles, requiring eight bits.

The direct enumeration is independent of the rank formula. It is implemented
in `succinct_profile_barrier.py` and included in `run_calibration.py`.

For the six-label DGY/HPR `Z_15` companion, both component Gram ranks are six.
Direct enumeration gives `3^6` and `5^6` component profiles, hence `15^6 =
11,390,625` joint profiles and a 24-bit lower bound. The materialized vector
dimension is 17, illustrating that the correct obstruction is true Gram rank,
not padded construction dimension.
