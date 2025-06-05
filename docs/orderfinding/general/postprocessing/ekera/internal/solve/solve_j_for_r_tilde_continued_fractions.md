## Function: <code>solve\_j\_for\_r\_tilde\_continued\_fractions(j, m, l)</code>
For j = j0(z) an optimal frequency, for m such that r < 2^m and l such that r^2 < 2^(m+l), and any z in [0, r), this function recovers r_tilde = r / d where d = gcd(r, z) by expanding the quotient j / 2^(m+l) in continued fractions and returning the last denominator < 2^((m+l)/2) as described in [[E24]](https://doi.org/10.1145/3655026).

By Lem. 4.1 in [[E24]](https://doi.org/10.1145/3655026), this function is guaranteed to return r_tilde provided that the requirements on the input parameters are met.

For further details, see Lem. 4.1, and Sect. 4 and App. B, of [[E24]](https://doi.org/10.1145/3655026).

## Import directive
```python
from quaspy.orderfinding.general.postprocessing.ekera.internal.solve import solve_j_for_r_tilde_continued_fractions
```

## Parent module
- [<code>solve</code>](README.md)

## Prototype
```python
def solve_j_for_r_tilde_continued_fractions(j,
                                            m,
                                            l)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| j | An optimal frequency j0(z), for m and l as passed to this function, and for any z in [0, r).<br><br>If some other frequency j on [0, 2^(m+l)) is passed to this function, it may or may not return r_tilde. |
| m | A positive integer m such that r < 2^m. |
| l | A positive integer l <= m, such that r^2 < 2^(m+l). |

## Return value
The last denominator < 2^((m+l)/2) in the continued fraction expansion of j / 2^(m+l). This denominator is guaranteed to be equal to r_tilde, provided that the requirements on the input parameters are met.

