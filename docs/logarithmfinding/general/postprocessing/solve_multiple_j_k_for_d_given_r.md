## Function: <code>solve\_multiple\_j\_k\_for\_d\_given\_r(j_k_list, m, sigma, l, g, x, r, ..)</code>
Attempts to compute the general discrete logarithm d = log_g x given a list of n frequency pairs [[j_1, k_1], ..., [j_n, k_n]] yielded by n independent runs of the quantum part of Shor's algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700) for computing general discrete logarithms, modified as in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084), and the order r of g, by using the lattice-based post-processing algorithm described in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084) (see Sect. 6) and [[E24t]](https://diva-portal.org/smash/get/diva2:1902626/FULLTEXT01.pdf) (see Sect. 5.5).

Note that this function does not implement meet-in-the-middle-techniques, although it is noted in [[E23p]](https://doi.org/10.48550/arXiv.2309.01754) that it is possible to use such techniques to speed up the post-processing also for general discrete logarithms.

## Import directive
```python
from quaspy.logarithmfinding.general.postprocessing import solve_multiple_j_k_for_d_given_r
```

## Parent module
- [<code>postprocessing</code>](README.md)

## Prototype
```python
def solve_multiple_j_k_for_d_given_r(j_k_list : list,
                                     m : int,
                                     sigma : int,
                                     l : int,
                                     g : CyclicGroupElement,
                                     x : CyclicGroupElement,
                                     r : int | gmpy2.mpz,
                                     tau : int = 0,
                                     delta : float = 0.99,
                                     precision : int | None = None,
                                     enumerate : bool | quaspy.logarithmfinding.general.postprocessing.EnumerationOptions = False,
                                     timeout : int | None | quaspy.utils.timeout.Timeout = None,
                                     verbose : bool = False)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| j_k_list | The n frequency pairs [[j_1, k_1], ..., [j_n, k_n]] where j_1, ..., j_n are integers on [0, 2^(m + sigma)) and k_1, ..., k_n are integers on [0, 2^l). |
| m | A positive integer m such that r < 2^m. |
| sigma | A non-negative integer sigma such that m + sigma is the length of the first control register in the quantum part of the algorithm. |
| l | A positive integer l ≈ m / s for s a tradeoff factor such that l is the length of the second control register in the quantum part of the algorithm. |
| g | The group element g. |
| x | The group element x = g^d. |
| r | The order r of g. |
| tau | A positive integer tau. Used to scale the basis for the lattice L^tau that is used in the post-processing. |
| delta | The parameter delta to use when delta-LLL-reducing the basis for the lattice L^tau used in the post-processing. Must be on the interval (1/4, 1]. A polynomial runtime in the dimension of the lattice is only guaranteed for delta < 1. |
| precision | The precision to use when computing the Gram–Schmidt projection factors as a part of delta-LLL-reducing the basis for the lattice L^tau used in the post-processing.<br><br>The precision may be set to None, as is the default, in which case the projection factors are represented as exact quotients. |
| enumerate | A flag that may be set to True to enumerate vectors in the lattice L^tau (until d is found or the specified timeout has elapsed), or to EnumerationOptions.CVP to consider only a closest vector in the lattice as returned by performing a limited enumeration, or to False to consider only the vector returned by Babai's nearest plane algorithm.<br><br>May also be set to EnumerationOptions.BOUNDED_BY_TAU in which case all vectors within distance R of the origin of the lattice L^tau are enumerated, where R depends on tau as R = sqrt(n + 1) * 2^(m + sigma - l + tau). |
| timeout | A timeout after which a TimeoutError will be raised and the computation aborted.<br><br>The timeout may be represented as an integer specifying the timeout in seconds, or as an instance of the Timeout class. May be set to None, as is the default, in which case no timeout is enforced. |
| verbose | A flag that may be set to True to print intermediary results and status updates when executing the post-processing algorithm. |

## Return value
The logarithm d, or None, if solving for d fails.

