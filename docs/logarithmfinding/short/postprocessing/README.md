## Module: <code>postprocessing</code>
A module for solving frequency pairs yielded by the quantum part of Ekerå–Håstad's algorithm [[EH17]](https://doi.org/10.1007/978-3-319-59879-6_20) for computing short discrete logarithms. This by using the post-processing algorithms in [[E20]](https://doi.org/10.1007/s10623-020-00783-2), [[E23p]](https://doi.org/10.48550/arXiv.2309.01754) and [[E24t]](https://diva-portal.org/smash/get/diva2:1902626/FULLTEXT01.pdf) (see Sect. 5.6).

## Import directive
```python
import quaspy.logarithmfinding.short.postprocessing
```

## Parent module
- [<code>short</code>](../README.md)

## Functions
- [<code>expected_u_for_j_k_d(j, k, m, l, d, ..)</code>](expected_u_for_j_k_d.md)

  Computes the vector u that is associated with a given frequency pair (j, k) and logarithm d, when using the post-processing algorithm from [[E23p]](https://doi.org/10.48550/arXiv.2309.01754) to post-process the output of the quantum part of Ekerå–Håstad's algorithm [[EH17]](https://doi.org/10.1007/978-3-319-59879-6_20).

- [<code>expected_u_for_multiple_j_k_d(j_k_list, m, l, d, ..)</code>](expected_u_for_multiple_j_k_d.md)

  Computes the vector u that is associated with a given list of n frequency pairs [[j_1, k_1], ..., [j_n, k_n]] and logarithm d, when using the post-processing algorithm from [[E20]](https://doi.org/10.1007/s10623-020-00783-2) or [[E23p]](https://doi.org/10.48550/arXiv.2309.01754) to post-process the output of the quantum part of Ekerå–Håstad's algorithm [[EH17]](https://doi.org/10.1007/978-3-319-59879-6_20).

- [<code>solve_j_k_for_d(j, k, m, l, g, x, tau, ..)</code>](solve_j_k_for_d.md)

  Attempts to compute the short discrete logarithm d = log_g x given a frequency pair (j, k) yielded by the quantum part of Ekerå–Håstad's algorithm [[EH17]](https://doi.org/10.1007/978-3-319-59879-6_20) by using the post-processing algorithms in [[E23p]](https://doi.org/10.48550/arXiv.2309.01754).

- [<code>solve_multiple_j_k_for_d(j_k_list, m, l, g, x, ..)</code>](solve_multiple_j_k_for_d.md)

  Attempts to compute the short discrete logarithm d = log_g x given a list of n frequency pairs [[j_1, k_1], ..., [j_n, k_n]] yielded by n independent runs of the quantum part of Ekerå–Håstad's algorithm [[EH17]](https://doi.org/10.1007/978-3-319-59879-6_20) by using the classical post-processing algorithm from [[E20]](https://doi.org/10.1007/s10623-020-00783-2) and [[E24t]](https://diva-portal.org/smash/get/diva2:1902626/FULLTEXT01.pdf) (see Sect. 5.6).

## Classes
- [<code>EnumerationOptions</code>](EnumerationOptions.md)

  An enumeration of options for solving a list of n frequency pairs [[j_1, k_1], ..., [j_n, k_n]] for d.
