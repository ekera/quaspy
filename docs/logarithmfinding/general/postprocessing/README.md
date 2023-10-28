## Module: <code>postprocessing</code>
A module for solving a frequency pair (j, k) yielded by the quantum part of Shor's algorithm for the general discrete logarithm problem for the logarithm d given the order r.

## Import directive
```python
import quaspy.logarithmfinding.general.postprocessing
```

## Parent module
- [<code>general</code>](../README.md)

## Functions
- [<code>solve_j_k_for_d_given_r(j, k, m, sigma, l, g, x, r, ..)</code>](solve_j_k_for_d_given_r.md)

  Attempts to compute the general discrete logarithm d given a frequency pair (j, k) yielded by the quantum part of Shor's algorithm as modified in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084), and the order r, by using the modified post-processing algorithm described in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084).

