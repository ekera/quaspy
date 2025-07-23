## Function: <code>solve\_j\_for\_r(j, m, l, g, ..)</code>
Attempts to compute the order r of g, or a positive integer multiple thereof, given a frequency j yielded by the quantum part of Shor's order-finding algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700), by using the post-processing algorithms described in detail in [[E24]](https://doi.org/10.1145/3655026).

The idea is to try to solve not only j, but also j ± 1, ..., j ± B, for r, with the aim of solving an optimal frequency j0(z) for r, for z the peak index on [0, r). Provided

- that j0(z) is solved for r,

- that d = gcd(r, z) is cm-smooth (for the definition of cm-smooth in [[E24]](https://doi.org/10.1145/3655026)),

- that l is selected as required by the solution method (see the documentation for parameter l below), and

- that the accept_multiple flag is set to False,

the order r will be found by this function.

Note that this function does not implement meet-in-the-middle-techniques, although it is noted in [[E23p]](https://doi.org/10.48550/arXiv.2309.01754) that it is possible to use such techniques to speed up the post-processing also for order finding.

## Import directive
```python
from quaspy.orderfinding.general.postprocessing.ekera import solve_j_for_r
```

## Parent module
- [<code>ekera</code>](README.md)

## Prototype
```python
def solve_j_for_r(j : int | gmpy2.mpz,
                  m : int,
                  l : int,
                  g : CyclicGroupElement,
                  c : int = 1,
                  B : int = 1000,
                  accept_multiple : bool = False,
                  method : SolutionMethods = SolutionMethods.LATTICE_BASED_SHORTEST_VECTOR,
                  timeout : int | None | quaspy.utils.timeout.Timeout = None,
                  verbose : bool = False,
                  opt_speculative : bool = True,
                  opt_isolate_peak : bool = True)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| j | The frequency j yielded by the quantum order-finding algorithm. |
| m | A positive integer m such that r < 2^m. |
| l | A positive integer l <= m, such that m + l is the length of the control register in the quantum order-finding algorithm.<br><br>If method is set to SolutionMethods.CONTINUED_FRACTIONS_BASED or SolutionMethods.LATTICE_BASED_SHORTEST_VECTOR, it is required that r^2 < 2^(m + l) or else r may not be found.<br><br>If method is set to SolutionMethods.LATTICE_BASED_ENUMERATE, it is possible to select l = m - Delta for some Delta in [0, m), at the expense of enumerating at most 6 * sqrt(3) * 2^Delta lattice vectors for each offset in j considered. |
| g | The group element g of order r. |
| c | A parameter c >= 1 that specifies the maximum size of the missing cm-smooth component d in r = d * r_tilde when solving j for r, for cm-smooth as defined in [[E24]](https://doi.org/10.1145/3655026).<br><br>As is explained in [[E24]](https://doi.org/10.1145/3655026), increasing c increases the success probability, at the expense of increasing the runtime. |
| B | A bound B >= 0 on the offset in j. If B > 0, this function tries to solve not only j, but also j ± 1, ..., j ± B, for r, or for a positive integer multiple of r. |
| accept_multiple | A flag that may be set to True to indicate that only a positive integer multiple of r is sought. If set to True, this function returns as soon as it finds r such that g^r = 1. |
| method | An enumeration entry from the SolutionMethods class that specifies the method to use to solve j for r. For further details, see the documentation for the SolutionMethods class. |
| timeout | A timeout after which a TimeoutError will be raised and the computation aborted.<br><br>The timeout may be represented as an integer specifying the timeout in seconds, or as an instance of the Timeout class. May be set to None, as is the default, in which case no timeout is enforced. |
| verbose | A flag that may be set to True to print intermediary results and status updates when executing the post-processing algorithm. |
| opt_speculative | A flag that may be set to True to indicate that Alg. 2 in [[E24]](https://doi.org/10.1145/3655026) should be used instead of Alg. 3 to find the missing cm-smooth component of r. In most cases, Alg. 2 is faster than Alg. 3, but in the worst case Alg. 2 is a lot slower than Alg. 3. For further details, see [[E24]](https://doi.org/10.1145/3655026). |
| opt_isolate_peak | A flag that may be set to True to indicate that all offsets in j up to B should not be tested exhaustively. Rather, the peak around the optimal frequency j_0(z) should be isolated: As soon as offsets to the left and right of j_0(z) are found for which the post-processing algorithm fails to produce r such that g^r = 1, this function returns the minimum r found such that g^r = 1. Note that this flag has no effect if the accept_multiple flag is set to True, as this function then returns as soon as it finds r such that g^r = 1. |

## Return value
If the accept_multiple flag is set to False, the order r of g is returned with probability >= P, for P as given by the lower bound in [[E24]](https://doi.org/10.1145/3655026) (provided that the opt_isolate_peak flag is set to False). Otherwise, None, or exceptionally a positive integer multiple of the order r, is returned. If the accept_multiple flag is set to True, some positive integer multiple of r is returned with probabilty >= P. Otherwise, None is returned.

