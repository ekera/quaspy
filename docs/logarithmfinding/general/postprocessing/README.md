## Module: <code>postprocessing</code>
A module for solving frequency pairs yielded by the quantum part of Shor's algorithm for the general discrete logarithm problem [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700), modified as in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084) and [[E24t]](https://diva-portal.org/smash/get/diva2:1902626/FULLTEXT01.pdf) (see Sect. 5.5), for the logarithm given the order.

## Import directive
```python
import quaspy.logarithmfinding.general.postprocessing
```

## Parent module
- [<code>general</code>](../README.md)

## Functions
- [<code>solve_j_k_for_d_given_r(j, k, m, sigma, l, g, x, r, ..)</code>](solve_j_k_for_d_given_r.md)

  Attempts to compute the general discrete logarithm d = log_g x given a frequency pair (j, k) yielded by the quantum part of Shor's algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700), modified as in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084), and the order r of g, by using the post-processing algorithm from [[E19p]](https://doi.org/10.48550/arXiv.1905.09084).

- [<code>solve_multiple_j_k_for_d_given_r(j_k_list, m, sigma, l, g, x, r, ..)</code>](solve_multiple_j_k_for_d_given_r.md)

  Attempts to compute the general discrete logarithm d = log_g x given a list of n frequency pairs [[j_1, k_1], ..., [j_n, k_n]] yielded by n independent runs of the quantum part of Shor's algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700) for computing general discrete logarithms, modified as in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084), and the order r of g, by using the lattice-based post-processing algorithm described in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084) (see Sect. 6) and [[E24t]](https://diva-portal.org/smash/get/diva2:1902626/FULLTEXT01.pdf) (see Sect. 5.5).

## Classes
- [<code>EnumerationOptions</code>](EnumerationOptions.md)

  An enumeration of options for solving a list of n frequency pairs [[j_1, k_1], ..., [j_n, k_n]] for d given r.
