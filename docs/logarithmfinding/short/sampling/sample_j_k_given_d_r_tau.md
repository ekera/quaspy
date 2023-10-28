## Function: <code>sample\_j\_k\_given\_d\_r\_tau(d, r, m, l, tau, ..)</code>
Samples a frequency pair (j, k) from the distribution induced by Ekerå–Håstad's quantum algorithm for finding a given short discrete logarithm d in a group of unknown order r.

To sample the distribution, j is first picked uniformly at random from [0, 2^(m + l)). A pivot is then selected uniformly at random from [0, 1), and the optimal frequency k0(j) for k computed.

For offsets 0, ±1, ±2, .., ±(B - 1) from k0(j) the probabilities 2^(m+l) * P(j, k = (k0(j) + offset) mod 2^l) of observing (j, k) are then subtracted from the pivot, and (j, k) returned as soon as pivot <= 0 provided that (j, k) is tau-good by Def. 1 in [[E23p]](https://doi.org/10.48550/arXiv.2309.01754). Otherwise, None is returned. The bound B is setup as a function of tau so that all tau-good pairs (j, k) are included in the search.

Note that it follows from the analysis in [[E20]](https://doi.org/10.1007/s10623-020-00783-2) and [[E23p]](https://doi.org/10.48550/arXiv.2309.01754) that j is distributed uniformly at random when r >= 2^(m + l) + (2^l - 1) * d. This function checks that this requirement is respected.

## Import directive
```python
from quaspy.logarithmfinding.short.sampling import sample_j_k_given_d_r_tau
```

## Parent module
- [<code>sampling</code>](README.md)

## Prototype
```python
def sample_j_k_given_d_r_tau(d,
                             r,
                             m,
                             l,
                             tau,
                             verbose = False,
                             extended_result = False)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| d | The discrete logarithm d in [1, r). |
| r | The order r. Used to check that r >= 2^(m + l) + (2^l - 1) * d. May be set to None in which case the check is not performed. |
| m | A positive integer m such that d < 2^m. |
| l | An non-negative integer l such that m + l is the length of the first control register in the quantum algorithm.<br><br>It is furthermore required that r >= 2^(m + l) + (2^l - 1) * d, as this function is based on the analysis in [[E20]](https://doi.org/10.1007/s10623-020-00783-2) that imposes this requirement so as to simplify the analysis. |
| tau | The parameter tau. An integer on (0, l]. |
| verbose | A flag that may be set to True to print intermediary results when sampling. |
| extended_result | A flag that may be set to True to not only return the frequency pair [j, k], but [[j, k], [k0(j), offset]]. |

## Return value
The frequency pair [j, k] sampled if the extended_result flag is set to False, or [[j, k], [k0(j), offset]] if the extended_result flag is set to True, or None if sampling failed because the upper bound B on the offset from the optimal frequency k0(j) was reached, or because (j, k) is not tau-good.

