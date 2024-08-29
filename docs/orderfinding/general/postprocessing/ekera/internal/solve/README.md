## Module: <code>solve</code>
A module for solving an optimal frequency j = j0(z) for r_tilde = r / d where d = gcd(r, z).

## Import directive
```python
import quaspy.orderfinding.general.postprocessing.ekera.internal.solve
```

## Parent module
- [<code>internal</code>](../README.md)

## Functions
- [<code>solve_j_for_r_tilde_continued_fractions(j, m, l)</code>](solve_j_for_r_tilde_continued_fractions.md)

  For j = j0(z) an optimal frequency, for m such that r < 2^m and l such that r^2 < 2^(m+l), and any z in [0, r), this function recovers r_tilde = r / d where d = gcd(r, z) by expanding the quotient j / 2^(m+l) in continued fractions and returning the last denominator < 2^((m+l)/2) as described in [[E24]](https://doi.org/10.1145/3655026).

- [<code>solve_j_for_r_tilde_lattice_enumerate(j, m, l, g, ..)</code>](solve_j_for_r_tilde_lattice_enumerate.md)

  For j = j0(z) an optimal frequency, for m such that r < 2^m and l = m - Delta for some Delta in [0, m), and z in [0, r), this function recovers r_tilde = r / d where d = gcd(r, z), provided that d is cm-smooth, by enumerating at most 6 * sqrt(3) * 2^Delta) vectors in a two-dimensional lattice L as described in [[E24]](https://doi.org/10.1145/3655026).

- [<code>solve_j_for_r_tilde_lattice_svp(j, m, l, ..)</code>](solve_j_for_r_tilde_lattice_svp.md)

  For j = j0(z) an optimal frequency, for m such that r < 2^m and l such that r^2 < 2^(m+l), and any z in [0, r), this function recovers r_tilde = r / d where d = gcd(r, z) by finding the shortest non-zero vector in a two-dimensional lattice L as described in [[E24]](https://doi.org/10.1145/3655026).

