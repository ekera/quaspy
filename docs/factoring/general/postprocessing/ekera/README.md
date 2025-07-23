## Module: <code>ekera</code>
A module for factoring N completely given the order r of an element g selected uniformly at random from the multiplicative group of the ring of integers modulo N, where g need not be explicitly specified. This by using the algorithm from [[E21b]](https://doi.org/10.1007/s11128-021-03069-1).

This module furthermore contains convenience functions for first solving the frequency j yielded by the quantum part of Shor's order-finding algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700) for a positive integer multiple r' of r, and for then solving r' for the complete factorization of N. This by using the classical post-processing algorithms from [[E24]](https://doi.org/10.1145/3655026) in the first step, and the algorithm from [[E21b]](https://doi.org/10.1007/s11128-021-03069-1) in the second step.

Finally, this module contains convenience functions for first solving a list of n frequencies [j_1, ..., j_n] yielded by n runs of the quantum part of Seifert's variation [[Seifert01]](https://doi.org/10.1007/3-540-45353-9_24) of Shor's order-finding algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700) for a positive integer multiple r' of r, and for then solving r' for the complete factorization of N. This by using the classical post-processing algorithms from [[E24t]](https://diva-portal.org/smash/get/diva2:1902626/FULLTEXT01.pdf) (see Sect. 5.4) and [[E21]](https://doi.org/10.1515/jmc-2020-0006), with supporting functions from [[E24]](https://doi.org/10.1145/3655026), in the first step, and by using the algorithm from [[E21b]](https://doi.org/10.1007/s11128-021-03069-1) in the second step.

Throughout this module, the algorithms are as described in [[E24t]](https://diva-portal.org/smash/get/diva2:1902626/FULLTEXT01.pdf), [[E24]](https://doi.org/10.1145/3655026), [[E21]](https://doi.org/10.1515/jmc-2020-0006) and [[E21b]](https://doi.org/10.1007/s11128-021-03069-1). The notation is also inherited from said works.

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

  Attempts to factor N completely, given the frequency j yielded by a run of the quantum part of Shor's order-finding algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700) when computing the order r of g modulo N, for g an element selected uniformly at random from the multiplicative group of the ring of integers modulo N.

- [<code>solve_j_for_factors_mod_N(j, m, l, g, N, ..)</code>](solve_j_for_factors_mod_N.md)

  Attempts to factor N completely, given the frequency j yielded by a run of the quantum part of Shor's order-finding algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700) when computing the order r of g modulo N, for g an element selected uniformly at random from the multiplicative group of the ring of integers modulo N.

- [<code>solve_multiple_j_for_factors(j_list, m, l, g, N, ..)</code>](solve_multiple_j_for_factors.md)

  Attempts to factor N completely, given a list of n frequencies [j_1, ..., j_n] yielded by n independent runs of the quantum part of Seifert's variation [[Seifert01]](https://doi.org/10.1007/3-540-45353-9_24) of Shor's order-finding algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700) when computing the order r of g modulo N, for g an element selected uniformly at random from the multiplicative group of the ring of integers modulo N.

- [<code>solve_multiple_j_for_factors_mod_N(j_list, m, l, g, N, ..)</code>](solve_multiple_j_for_factors_mod_N.md)

  Attempts to factor N completely, given a list of n frequencies [j_1, ..., j_n] yielded by n independent runs of the quantum part of Seifert's variation [[Seifert01]](https://doi.org/10.1007/3-540-45353-9_24) of Shor's order-finding algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700) when computing the order r of g modulo N, for g an element selected uniformly at random from the multiplicative group of the ring of integers modulo N.

- [<code>solve_r_for_factors(r, N, ..)</code>](solve_r_for_factors.md)

  Attempts to factor N completely given the order r of an element g selected uniformly at random from the multiplicative group of the ring of integers modulo N, where g need not be explicitly specified. This by using the algorithm from [[E21b]](https://doi.org/10.1007/s11128-021-03069-1).

## Classes
- [<code>IncompleteFactorizationException</code>](IncompleteFactorizationException.md)

  An exception that is raised to signal an incomplete factorization.
- [<code>OptProcessCompositeFactors</code>](OptProcessCompositeFactors.md)

  An enumeration of optimization options for how the solver is to select x, and to process composite factors.
