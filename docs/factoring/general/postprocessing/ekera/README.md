## Module: <code>ekera</code>
A module for factoring N completely given the order r of an element g selected uniformly at random from the multiplicative group of the ring of integers modulo N, where g need not be explicitly specified. This by using the algorithm in [[E21b]](https://doi.org/10.1007/s11128-021-03069-1).

This module furthermore contains convenience functions for first solving the frequency j yielded by Shor's order-finding algorithm for a positive integer multiple r' of r, and for then solving r' for the complete factorization of N. This by using the algorithms in [[E22p]](https://doi.org/10.48550/arXiv.2201.07791) in the first step, and the algorithm in [[E21b]](https://doi.org/10.1007/s11128-021-03069-1) in the second step.

## Import directive
```python
import quaspy.factoring.general.postprocessing.ekera
```

## Parent module
- [<code>postprocessing</code>](../README.md)

## Submodules
- [<code>internal</code>](internal/README.md)

  A module for internal classes and functions used by the functions in the parent module for factoring N completely.

## Functions
- [<code>solve_j_for_factors(j, m, l, g, N, ..)</code>](solve_j_for_factors.md)

  Attempts to factor N completely, given the frequency j yielded by the quantum order-finding algorithm when computing the order r of g modulo N, for g an element selected uniformly at random from the multiplicative group of the ring of integers modulo N.

- [<code>solve_j_for_factors_mod_N(j, m, l, g, N, ..)</code>](solve_j_for_factors_mod_N.md)

  Attempts to factor N completely, given the frequency j yielded by the quantum order-finding algorithm when computing the order r of g modulo N, for g an element selected uniformly at random from the multiplicative group of the ring of integers modulo N.

- [<code>solve_r_for_factors(r, N, ..)</code>](solve_r_for_factors.md)

  Attempts to factor N completely given the order r of an element g selected uniformly at random from the multiplicative group of the ring of integers modulo N, where g need not be explicitly specified.

## Classes
- [<code>IncompleteFactorizationException</code>](IncompleteFactorizationException.md)

  An exception that is raised to signal an incomplete factorization.
- [<code>OptProcessCompositeFactors</code>](OptProcessCompositeFactors.md)

  An enumeration of optimization options for how the solver is to select x, and to process composite factors.
