## Function: <code>solve\_j\_for\_factors\_mod\_N(j, m, l, g, N, ..)</code>
Attempts to factor N completely, given the frequency j yielded by the quantum order-finding algorithm when computing the order r of g modulo N, for g an element selected uniformly at random from the multiplicative group of the ring of integers modulo N.

This by using the algorithm in [[E21b]](https://doi.org/10.1007/s11128-021-03069-1) to factor N given r, or a positive multiple of r, and the post-processing algorithm in [[E22p]](https://doi.org/10.48550/arXiv.2201.07791) to find r, or a positive multiple of r, given j.

> This convenience function simply calls solve_j_for_factors() with g setup by calling IntegerModRingMulSubgroupElement(g, N).

> In turn, the solve_j_for_factors() convenience function simply calls solve_j_for_r(), and then solve_r_for_factors(), passing along r. To access all options of these functions, call them manually in sequence instead.

## Import directive
```python
from quaspy.factoring.general.postprocessing.ekera import solve_j_for_factors_mod_N
```

## Parent module
- [<code>ekera</code>](README.md)

## Prototype
```python
def solve_j_for_factors_mod_N(j : int,
                              m : int,
                              l : int,
                              g : int,
                              N : int,
                              c_solve : int = 1,
                              c_factor : int = 1,
                              B : int = 1000,
                              k = None,
                              timeout = None,
                              accept_multiple = True,
                              method = SolutionMethods.LATTICE_BASED_SHORTEST_VECTOR,
                              verbose = False,
                              opt_speculative = True)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| j | The frequency j yielded by the quantum order-finding algorithm. |
| m | A positive integer m such that r < 2^m. |
| l | A positive integer l <= m, such that m+l is the length of the control register in the quantum order-finding algorithm.<br><br>If method is set to SolutionMethods.CONTINUED_FRACTIONS_BASED or SolutionMethods.LATTICE_BASED_SHORTEST_VECTOR, it is required that r^2 < 2^(m+l) or else r may not be found.<br><br>If method is set to SolutionMethods.LATTICE_BASED_ENUMERATE, it is possible to select l = m - Delta for some Delta in [0, m), at the expense of enumerating at most 6 * sqrt(3) * 2^Delta lattice vectors for each offset in j considered. |
| g | The group element g of order r. |
| N | The integer N. |
| c_solve | A parameter c_solve >= 1 that specifies the maximum size of the missing smooth component d in r = d * r_tilde when solving j for r, or a multiple of r. As is explained in [[E22p]](https://doi.org/10.48550/arXiv.2201.07791), increasing c increases the success probability, at the expense of increasing the runtime. |
| c_factor | A parameter c_factor >= 1 that specifies the maximum size of the missing smooth component in lambda'(N) when solving r for the complete factorization of N. As is explained in [[E21b]](https://doi.org/10.1007/s11128-021-03069-1), increasing c increases the success probability, at the expense of increasing the runtime. |
| B | A bound B >= 0 on the offset in j. If B > 0, the solve_j_for_r() function called by this function tries to solve not only j, but also j ± 1, .., j ± B, for r, or for a positive integer multiple of r. |
| k | The maximum number of iterations to perform when factoring N given the order r, or a multiple of r. Defaults to None.<br><br>If k is set to None, as many iterations as are necessary to completely factor N will be performed. If k is explicitly specified, and the complete factorization of N has not been found after k iterations, an exception of type IncompleteFactorizationException will be raised. |
| timeout | A timeout in seconds. Defaults to None. If the timeout is set to None, as much time as is necessary to completely factor N given the order r, or a multiple of r, will be used. If a timeout is explicitly specified, and the complete factorization of N has not been found when the timeout elapses, an exception of type IncompleteFactorizationException will be raised. |
| accept_multiple | A flag that may be set to True to indicate that only a positive integer multiple of r is sought. If set to True, the solve_j_for_r() function called by this function returns as soon as it finds r such that g^r = 1. |
| method | An enumeration entry from the SolutionMethods class that specifies the method to use to solve j for r. For further details, see the documentation for the SolutionMethods class. |
| verbose | A flag that may be set to True to print intermediary results and status updates when executing the post-processing algorithm. |
| opt_speculative | A flag that may be set to True to indicate that Algorithm 2 in [[E22p]](https://doi.org/10.48550/arXiv.2201.07791) should be used instead of Algorithm 3 to find the missing cm-smooth component of r. In most cases, Algorithm 2 is faster than Algorithm 3, but in the worst case Algorithm 2 is a lot slower than Algorithm 3. For further details, see [[E22p]](https://doi.org/10.48550/arXiv.2201.07791). |

## Return value
The set of all distinct prime factors that divide N, or None, if a positive multiple of r could not be found given j.

