## Function: <code>sample\_j\_k\_given\_d\_r(d, r, m, l, ..)</code>
Samples a frequency pair (j, k) from the distribution induced by the quantum part of Ekerå–Håstad's algorithm [[EH17]](https://doi.org/10.1007/978-3-319-59879-6_20) for finding a short discrete logarithm d in a group of unknown order r.

The sampling procedure is described in Sect. 5.6.5 of [[E24t]](https://diva-portal.org/smash/get/diva2:1902626/FULLTEXT01.pdf).

To sample the distribution, j is first picked uniformly at random from [0, 2^(m + l)). A pivot is then selected uniformly at random from [0, 1), and the optimal frequency k0(j) for k computed.

For offsets 0, ±1, ±2, ..., ±(B - 1) from k0(j) the probabilities 2^(m+l) * P(j, k = (k0(j) + offset) mod 2^l) of observing (j, k) are then subtracted from the pivot, and (j, k) returned as soon as pivot <= 0.

Note that it follows from the analysis in [[E20]](https://doi.org/10.1007/s10623-020-00783-2) and [[E23p]](https://doi.org/10.48550/arXiv.2309.01754) that j is distributed uniformly at random when r >= 2^(m + l) + (2^l - 1) * d. This function checks that this requirement is respected provided that r is included in the function call.

## Import directive
```python
from quaspy.logarithmfinding.short.sampling import sample_j_k_given_d_r
```

## Parent module
- [<code>sampling</code>](README.md)

## Prototype
```python
def sample_j_k_given_d_r(d : int | gmpy2.mpz,
                         r : int | gmpy2.mpz | None,
                         m : int,
                         l : int,
                         B : int = 100000,
                         timeout : int | None | quaspy.utils.timeout.Timeout = None,
                         verbose : bool = False,
                         extended_result : bool = False)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| d | The discrete logarithm d in [1, r). |
| r | The order r. Used to check that r >= 2^(m + l) + (2^l - 1) * d. May be set to None in which case the check is not performed. |
| m | A positive integer m such that d < 2^m. |
| l | A positive integer l such that the control registers in the quantum part of the algorithm are of lengths m + l and l qubits, respectively.<br><br>It is required that r >= 2^(m + l) + (2^l - 1) * d as this function is based on the analysis in [[E20]](https://doi.org/10.1007/s10623-020-00783-2) that imposes this requirement so as to simplify the analysis. |
| B | A parameter B that upper-bounds the offset from the optimal frequency k0(j) as explained above. |
| timeout | A timeout after which a TimeoutError will be raised and the sampling procedure aborted.<br><br>The timeout may be represented as an integer specifying the timeout in seconds, or as an instance of the Timeout class. May be set to None, as is the default, in which case no timeout is enforced. |
| verbose | A flag that may be set to True to print intermediary results when sampling. |
| extended_result | A flag that may be set to True to not only return the frequency pair (j, k) as a list [j, k], but [[j, k], [k0(j), offset]]. |

## Return value
The frequency pair [j, k] sampled if the extended_result flag is set to False, or [[j, k], [k0(j), offset]] if the extended_result flag is set to True, or None if sampling failed because the upper bound B on the offset from the optimal frequency k0(j) was reached.

