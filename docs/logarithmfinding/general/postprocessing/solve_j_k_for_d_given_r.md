## Function: <code>solve\_j\_k\_for\_d\_given\_r(j, k, m, sigma, l, g, x, r, ..)</code>
Attempts to compute the general discrete logarithm d = log_g x given a frequency pair (j, k) yielded by the quantum part of Shor's algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700), modified as in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084), and the order r of g, by using the post-processing algorithm from [[E19p]](https://doi.org/10.48550/arXiv.1905.09084).

More specifically, the post-processing algorithm is in Sect. 6 of [[E19p]](https://doi.org/10.48550/arXiv.1905.09084).

Note that this function does not implement meet-in-the-middle-techniques, although it is noted in [[E23p]](https://doi.org/10.48550/arXiv.2309.01754) that it is possible to use such techniques to speed up the post-processing also for general discrete logarithms.

## Import directive
```python
from quaspy.logarithmfinding.general.postprocessing import solve_j_k_for_d_given_r
```

## Parent module
- [<code>postprocessing</code>](README.md)

## Prototype
```python
def solve_j_k_for_d_given_r(j : int | gmpy2.mpz,
                            k : int | gmpy2.mpz,
                            m : int,
                            sigma : int,
                            l : int,
                            g : CyclicGroupElement,
                            x : CyclicGroupElement,
                            r : int | gmpy2.mpz,
                            B_ETA : int = 1000,
                            B_T : int = 1000,
                            timeout : int | None | quaspy.utils.timeout.Timeout = None,
                            verbose : bool = False)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| j | The frequency j. An integer on [0, 2^(m + sigma)). |
| k | The frequency k. An integer on [0, 2^l). |
| m | A positive integer m such that r < 2^m. |
| sigma | A non-negative integer sigma such that m + sigma is the length of the first control register in the quantum part of the algorithm. |
| l | A positive integer l such that l is the length of the second control register in the quantum part of the algorithm. |
| g | The group element g. |
| x | The group element x = g^d. |
| r | The order r of g. |
| B_ETA | An upper bound on the search space in eta. |
| B_T | An upper bound on the search space in t. |
| timeout | A timeout after which a TimeoutError will be raised and the computation aborted.<br><br>The timeout may be represented as an integer specifying the timeout in seconds, or as an instance of the Timeout class. May be set to None, as is the default, in which case no timeout is enforced.<br><br> |
| verbose | A flag that may be set to True to print intermediary results and status updates when executing the post-processing algorithm. |

## Return value
The logarithm d, or None, if solving for d fails.

