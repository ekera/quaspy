## Function: <code>sample\_j\_k\_given\_d\_r\_heuristic(d, r, m, sigma, l, ..)</code>
Samples a frequency pair (j, k) heuristically from the distribution induced by Shor's quantum algorithm for finding a given discrete logarithm d in a group of known order r, or from the distribution induced by Ekerå–Håstad's and Ekerå's quantum algorithms, depending on how parameters are selected.

The sampling procedure is described in Sect. 5 of [[E19p]](https://doi.org/10.48550/arXiv.1905.09084).

## Import directive
```python
from quaspy.logarithmfinding.sampling import sample_j_k_given_d_r_heuristic
```

## Parent module
- [<code>sampling</code>](README.md)

## Prototype
```python
def sample_j_k_given_d_r_heuristic(d,
                                   r,
                                   m,
                                   sigma,
                                   l,
                                   B_DELTA = 1000,
                                   B_ETA = 10000,
                                   verbose = False,
                                   extended_result = False)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| d | The discrete logarithm d in [1, r). |
| r | The order r. |
| m | A positive integer m such that r < 2^m when sampling from the distribution induced by Shor's algorithm or Ekerå's algorithm, and such that d < 2^m when sampling from the distribution induced by Ekerå-Håstad's algorithm. |
| sigma | An non-negative integer sigma such that m + sigma is the length of the first control register in the quantum algorithm. |
| l | A positive integer l such that l is the length of the second control register in the quantum algorithm. |
| B_DELTA | A parameter that upper-bounds the offset from the optimal frequency k0(j) when sampling k given j and eta. |
| B_ETA | A parameter that upper-bounds eta when sampling j and eta. |
| verbose | A flag that may be set to True to print intermediary results when sampling. |
| extended_result | A flag that may be set to True to not only return the frequency pair [j, k], but [[j, k], [k0(j), offset, eta]]. |

## Return value
The frequency pair [j, k] sampled if the extended_result flag is set to False, or [[j, k], [k0(j), offset, eta]] if the extended_result flag is set to True, or None if sampling failed because the upper bound on the offset from the optimal frequency k0(j) or on eta were reached.

