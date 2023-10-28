## Class: <code>SolutionMethods</code>
An enumeration of methods for solving j for r.

As is explained in [[E22p]](https://doi.org/10.48550/arXiv.2201.07791), these methods are all designed to solve an optimal frequency j = j0(z), where z in [0, r), for r_tilde = r / d, where d = gcd(r, z) and d is cm-smooth. In what follows below, it is assumed that these requirements are met.

There are three solution methods:

1. CONTINUED_FRACTIONS_BASED

Expand j / 2^(m+l) in continued fractions to find z / r, and hence r_tilde = r / gcd(r, z), as originally proposed in [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700), but with slightly smaller m+l, as decribed in [[E22p]](https://doi.org/10.48550/arXiv.2201.07791).

By Lemma 6 in [[E22p]](https://doi.org/10.48550/arXiv.2201.07791), the last convergent p / q with denominator q < 2^((m+l)/2) in the continued fractions expansion of j / 2^(m+l) is equal to z / r, provided that r < 2^m and r^2 < 2^(m+l).

For further details, see Lemma 6, and Sect. 4 and App. B, of [[E22p]](https://doi.org/10.48550/arXiv.2201.07791).

2. LATTICE_BASED_SHORTEST_VECTOR

Use Lagrange's lattice basis reduction algorithm to find the shortest non-zero vector in the lattice L spanned by (j, 1/2) and (2^(m+l), 0).

By Lemma 7 in [[E22p]](https://doi.org/10.48550/arXiv.2201.07791), provided that r < 2^m and r^2 < 2^(m+l), the second component of the shortest non-zero vector in L has r_tilde / 2 in its second component, up to sign of course.

For further details, see Lemma 7, and Sect. 4 and App. C, of [[E22p]](https://doi.org/10.48550/arXiv.2201.07791).

3. LATTICE_BASED_ENUMERATE

Use Lagrange's lattice basis reduction algorithm to find a reduced basis for the lattice L spanned by (j, 1/2) and (2^(m+l), 0), and enumerate all vectors within a ball of radius 2^(m-1/2) in L centered at the origin to find u = (rj - 2^(m+l) z, r / 2) / d, and hence r_tilde, as the second component is r_tilde / 2.

By Lemma 8 in [[E22p]](https://doi.org/10.48550/arXiv.2201.07791), provided that r < 2^m and l = m - Delta, at most 6 * sqrt(3) * 2^Delta vectors must be enumerated in L to find u and hence r_tilde, so if Delta is small then this method is efficient.

In practice, as mentioned in [[E22p]](https://doi.org/10.48550/arXiv.2201.07791), the leading constant in the above bound is not tight, and the enumeration can be optimized. Some of these optimizations are implemented here so the enumeration typically considers fewer vectors than the bound indicates.

For further details, see Lemma 8, and Sect. 4 and App. C, of [[E22p]](https://doi.org/10.48550/arXiv.2201.07791).

## Import directive
```python
from quaspy.orderfinding.general.postprocessing.ekera import SolutionMethods
```

## Parent module
- [<code>ekera</code>](README.md)

## Members
- <code>CONTINUED_FRACTIONS_BASED</code>
- <code>LATTICE_BASED_SHORTEST_VECTOR</code>
- <code>LATTICE_BASED_ENUMERATE</code>

