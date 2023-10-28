## Module: <code>postprocessing</code>
A module for solving a frequency pair (j, k) yielded by the quantum part of Ekerå–Håstad's quantum algorithm for the order r. This by using the post-processing algorithms in [[E20]](https://doi.org/10.1007/s10623-020-00783-2) and [[E23p]](https://doi.org/10.48550/arXiv.2309.01754).

This implementation does not currently support tradeoffs, as it uses Lagrange's lattice basis reduction algorithm.

## Import directive
```python
import quaspy.logarithmfinding.short.postprocessing
```

## Parent module
- [<code>short</code>](../README.md)

## Functions
- [<code>expected_u_for_j_k_d(j, k, m, l, d, tau)</code>](expected_u_for_j_k_d.md)

  Computes the vector u that we seek to find for a given frequency pair (j, k) and logarithm d, when using the post-processing algorithms described in [[E20]](https://doi.org/10.1007/s10623-020-00783-2) and [[E23p]](https://doi.org/10.48550/arXiv.2309.01754).

- [<code>solve_j_k_for_d(j, k, m, l, g, x, tau, ..)</code>](solve_j_k_for_d.md)

  Attempts to compute the short discrete logarithm d given a frequency pair (j, k) yielded by the quantum part of Ekerå–Håstad's algorithm, by using the post-processing algorithms described in [[E23p]](https://doi.org/10.48550/arXiv.2309.01754).

