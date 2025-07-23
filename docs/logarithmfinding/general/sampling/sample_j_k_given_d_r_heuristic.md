## Function: <code>sample\_j\_k\_given\_d\_r\_heuristic(d, r, m, sigma, l, ..)</code>
Samples a frequency pair (j, k) heuristically from the distribution induced by the quantum part of Shor's algorithm for finding a discrete logarithm d in a group of known order r, or from the distribution induced by the quantum part of Ekerå's variation thereof, depending on how parameters are selected.

The sampling procedure is described in Sect. 5 of [[E19p]](https://doi.org/10.48550/arXiv.1905.09084).

## Import directive
```python
from quaspy.logarithmfinding.general.sampling import sample_j_k_given_d_r_heuristic
```

## Parent module
- [<code>sampling</code>](README.md)

## Prototype
```python
def sample_j_k_given_d_r_heuristic(d : int | gmpy2.mpz,
                                   r : int | gmpy2.mpz,
                                   m : int,
                                   sigma : int,
                                   l : int,
                                   B_DELTA : int = 1000,
                                   B_ETA : int = 10000,
                                   integration_steps : int = 128,
                                   timeout : int | None | quaspy.utils.timeout.Timeout = None,
                                   verbose : bool = False,
                                   extended_result : bool = False)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| d | The discrete logarithm d in [1, r). |
| r | The order r. |
| m | A positive integer m such that r < 2^m. |
| sigma | A non-negative integer sigma such that m + sigma is the length of the first control register in the quantum part of the algorithm. |
| l | A positive integer l such that l is the length of the second control register in the quantum part of the algorithm. |
| B_DELTA | A parameter that upper-bounds the offset from the optimal frequency k0(j) when sampling k given j and eta. |
| B_ETA | A parameter that upper-bounds eta when sampling j and eta. |
| integration_steps | The number of steps to perform when integrating the probability distribution. |
| timeout | A timeout after which a TimeoutError will be raised and the sampling procedure aborted.<br><br>The timeout may be represented as an integer specifying the timeout in seconds, or as an instance of the Timeout class. May be set to None, as is the default, in which case no timeout is enforced. |
| verbose | A flag that may be set to True to print intermediary results when sampling. |
| extended_result | A flag that may be set to True to not only return the frequency pair (j, k) as a list [j, k], but [[j, k], [k0(j), offset, eta]]. |

## Return value
The frequency pair [j, k] sampled if the extended_result flag is set to False, or [[j, k], [k0(j), offset, eta]] if the extended_result flag is set to True, or None if sampling failed because the upper bound on the offset from the optimal frequency k0(j) or on eta were reached.

