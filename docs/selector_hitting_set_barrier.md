# Coefficient-menu barrier for the HPR selector

Status: general counting bound plus exact small-prime calibration, 1 August
2026. This applies to coefficient-extraction variants of HPR Lemma 3.7.

## Question

HPR retain a chosen nonzero monomial `(alpha,beta)` of

`X^(g1 g2)-X^((g1+1)g2)-X^(g1(g2+1))+X^((g1+1)(g2+1))`.

Their implementation stores `beta` as an element of `F_p`, with an
`O(log p)`-bit budget. Could one choose from a much smaller fixed menu of
exponents and store only the menu index?

## Counting barrier

For a fixed exponent `beta`, its coefficient can be nonzero only if `beta`
equals one of

`g1 g2`, `(g1+1)g2`, `g1(g2+1)`, `(g1+1)(g2+1)` modulo `p`.

Each bilinear equation covers `O(p)` of the `p^2` base pairs. Thus a fixed
exponent occurs in the support of only `O(p)` selector polynomials. Any menu
that intersects the support for all `p^2` base pairs must contain `Omega(p)`
exponents. Selecting a guaranteed nonzero coefficient from such a menu still
requires `Omega(log p)` persistent bits.

This conclusion is independent of the lexicographic convention. It applies to
every coefficient-extraction rule that first fixes a universal exponent menu
and then records which menu element is used.

## Exact finite calibration

`selector_hitting_set.py` constructs every support and exhausts the minimum
menu. The first values are:

| Prime | Minimum menu size | One witness |
|---:|---:|---|
| 3 | 2 | `{0,1}` |
| 5 | 3 | `{0,1,3}` |
| 7 | 4 | `{0,2,3,5}` |
| 11 | 7 | `{0,1,2,3,5,7,9}` |
| 13 | 7 | `{1,2,4,6,8,10,12}` |

The witnesses are not canonical; only the minimum cardinalities matter.

## Consequence

The projected selector avoids `beta_k` at its output prime by replacing
coefficient extraction with the exact mixed derivative in the same
characteristic. The counting bound explains why merely optimizing the
remaining coefficient choices cannot remove their `log(m/p_k)` state. A full
recursive improvement must batch across characteristics or change the
reconstruction interface, rather than select monomials more cleverly.
