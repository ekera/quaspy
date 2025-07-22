## Function: <code>solve\_j\_for\_factors(j, m, l, g, N, ..)</code>
Attempts to factor N completely, given the frequency j yielded by a run of the quantum part of Shor's order-finding algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700) when computing the order r of g modulo N, for g an element selected uniformly at random from the multiplicative group of the ring of integers modulo N.

This by using the algorithm from [[E21b]](https://doi.org/10.1007/s11128-021-03069-1) to factor N given r, or a positive multiple of r, and the post-processing algorithm from [[E24]](https://doi.org/10.1145/3655026) to find r, or a positive multiple of r, given j.

Throughout this function, the algorithms are as described in [[E24]](https://doi.org/10.1145/3655026) and [[E21b]](https://doi.org/10.1007/s11128-021-03069-1). The notation is also inherited from said works.

> This convenience function simply calls solve_j_for_r(), and then solve_r_for_factors(), passing along r. To access all options of these functions, call them manually in sequence instead.

## Import directive
```python
from quaspy.factoring.general.postprocessing.ekera import solve_j_for_factors
```

## Parent module
- [<code>ekera</code>](README.md)

## Prototype
```python
def solve_j_for_factors(j : int | gmpy2.mpz,
                        m : int,
                        l : int,
                        g : CyclicGroupElement,
                        N : int | gmpy2.mpz,
                        c_solve : int = 1,
                        c_factor : int = 1,
                        B : int = 1000,
                        k : int | None = None,
                        accept_multiple : bool = True,
                        method : SolutionMethods = SolutionMethods.LATTICE_BASED_SHORTEST_VECTOR,
                        timeout : int | None | quaspy.utils.timeout.Timeout = None,
                        verbose : bool = False,
                        opt_speculative : bool = True)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| j | The frequency j yielded by the quantum part of the order-finding algorithm. |
| m | A positive integer m such that r < 2^m. |
| l | A positive integer l <= m such that m + l is the length of the control register in the quantum part of the order-finding algorithm.<br><br>If method is set to SolutionMethods.CONTINUED_FRACTIONS_BASED or SolutionMethods.LATTICE_BASED_SHORTEST_VECTOR, it is required that r^2 < 2^(m + l) or else r may not be found.<br><br>If method is set to SolutionMethods.LATTICE_BASED_ENUMERATE, it is possible to select l = m - Delta for some Delta in [0, m), at the expense of enumerating at most 6 * sqrt(3) * 2^Delta lattice vectors for each offset in j considered. |
| g | The group element g of order r. |
| N | The integer N. |
| c_solve | A parameter c_solve >= 1 that specifies the maximum size of the missing smooth component d in r = d * r_tilde when solving j for r, or a multiple of r. As is explained in [[E24]](https://doi.org/10.1145/3655026), increasing c increases the success probability, at the expense of increasing the runtime. |
| c_factor | A parameter c_factor >= 1 that specifies the maximum size of the missing smooth component in lambda'(N) when solving r for the complete factorization of N. As is explained in [[E21b]](https://doi.org/10.1007/s11128-021-03069-1), increasing c increases the success probability, at the expense of increasing the runtime. |
| B | A bound B >= 0 on the offset in j. If B > 0, the solve_j_for_r() function called by this function tries to solve not only j, but also j ± 1, ..., j ± B, for r, or for a positive integer multiple of r. |
| k | The maximum number of iterations to perform when factoring N given the order r, or a multiple of r. Defaults to None.<br><br>If k is set to None, as many iterations as are necessary to completely factor N will be performed. If k is explicitly specified, and the complete factorization of N has not been found after k iterations, an exception of type IncompleteFactorizationException will be raised. |
| accept_multiple | A flag that may be set to True to indicate that only a positive integer multiple of r is sought. If set to True, the solve_j_for_r() function called by this function returns as soon as it finds r such that g^r = 1. |
| method | An enumeration entry from the SolutionMethods class that specifies the method to use to solve j for r. For further details, see the documentation for the SolutionMethods class. |
| timeout | A timeout after which an IncompleteFactorizationException or TimeoutError will be raised and the computation aborted. More specifically, if the process of factoring N given r, or a multiple of r, has been initiated, an IncompleteFactorizationException will be raised, otherwise a TimeoutError will be raised, when the timeout elapses.<br><br>The timeout may be represented as an integer specifying the timeout in seconds, or as an instance of the Timeout class. May be set to None, as is the default, in which case no timeout is enforced. |
| verbose | A flag that may be set to True to print intermediary results and status updates when executing the post-processing algorithm. |
| opt_speculative | A flag that may be set to True to indicate that Alg. 2 in [[E24]](https://doi.org/10.1145/3655026) should be used instead of Alg. 3 to find the missing cm-smooth component of r. In most cases, Alg. 2 is faster than Alg. 3, but in the worst case Alg. 2 is a lot slower than Alg. 3. For further details, see [[E24]](https://doi.org/10.1145/3655026). |

## Return value
The set of all distinct prime factors that divide N, or None, if a positive multiple of r could not be found given j.

