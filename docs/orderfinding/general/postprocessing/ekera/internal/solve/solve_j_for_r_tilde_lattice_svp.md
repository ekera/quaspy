## Function: <code>solve\_j\_for\_r\_tilde\_lattice\_svp(j, m, l, ..)</code>
For j = j0(z) an optimal frequency, for m such that r < 2^m and l such that r^2 < 2^(m+l), and any z in [0, r), this function recovers r_tilde = r / d where d = gcd(r, z) by finding the shortest non-zero vector in a two-dimensional lattice L as described in [[E24]](https://doi.org/10.1145/3655026).

More specifically, this function uses Lagrange's lattice basis reduction algorithm to find the shortest non-zero vector

u = (rj - 2^(m+l) z, r / 2) / d

in the lattice L spanned by (j, 1/2) and (2^(m+l), 0), and hence r_tilde, as the second component is r_tilde / 2. This function return r_tilde.

By Lem. 4.2 in [[E24]](https://doi.org/10.1145/3655026), provided that r < 2^m and r^2 < 2^(m+l), the second component of the shortest non-zero vector in L has r_tilde / 2 in its second component, up to sign of course.

For further details, see Lem. 4.2, and Sect. 4 and App. C, of [[E24]](https://doi.org/10.1145/3655026).

## Import directive
```python
from quaspy.orderfinding.general.postprocessing.ekera.internal.solve import solve_j_for_r_tilde_lattice_svp
```

## Parent module
- [<code>solve</code>](README.md)

## Prototype
```python
def solve_j_for_r_tilde_lattice_svp(j,
                                    m,
                                    l,
                                    multiples = None)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| j | An optimal frequency j0(z), for m and l as passed to this function, and for any z in [0, r).<br><br>If some other frequency j on [0, 2^(m+l)) is passed to this function, it will also return an integer, that may or may not be equal to r_tilde. |
| m | A positive integer m such that r < 2^m. |
| l | A positive integer l <= m, such that r^2 < 2^(m+l). |
| multiples | Row multiples to use as a starting point when computing the Lagrange-reduced basis for the slightly scaled basis A = [[j, 1], [2^(m+l+1), 0]], see the lagrange() function for further details, or None, as is the default, if no such multiples are available.<br><br>The idea is that when trying to solve not only j but also j ± 1, .., j ± B for r_tilde, in the hope that this will lead to an optimal frequency j0(z) being solved for r_tilde, the row multiples that yield a Lagrange-reduced basis for adjecent offsets in j are likely to be close.<br><br>For this reason, this function accepts row multiples as input, and it furthermore returns as output the row multiples that yield a Lagrange-reduced basis for the value of j input, so that these can be fed back to this function when processing j + 1 or j - 1, recursively. |

## Return value
The tuple [r_tilde_candidate, multiple'], where multiple' are the row multiples that yield a reduced basis for the value of j input, and r_tilde_candidate is equal to r_tilde provided that the requirements on the input parameters are met.

