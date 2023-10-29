## Function: <code>solve\_j\_k\_for\_d\_given\_r(j, k, m, sigma, l, g, x, r, ..)</code>
Attempts to compute the general discrete logarithm d given a frequency pair (j, k) yielded by the quantum part of Shor's algorithm as modified in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084), and the order r, by using the modified post-processing algorithm described in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084).

The modified post-processing algorithm is described in Sect. 6 of [[E19p]](https://doi.org/10.48550/arXiv.1905.09084):

Note that this function does not implement meet-in-the-middle-techniques, although it is noted in [[E23p]](https://doi.org/10.48550/arXiv.2309.01754) that it is possible to use such techniques to speed up the post-processing also for general discrete logarithms.

## Import directive
```python
from quaspy.logarithmfinding.general.postprocessing import solve_j_k_for_d_given_r
```

## Parent module
- [<code>postprocessing</code>](README.md)

## Prototype
```python
def solve_j_k_for_d_given_r(j,
                            k,
                            m,
                            sigma,
                            l,
                            g : CyclicGroupElement,
                            x : CyclicGroupElement,
                            r,
                            B_ETA = 1000,
                            B_T = 1000,
                            verbose = False)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| j | The frequency j. |
| k | The frequency k. |
| m | A positive integer m such that d < 2^m. |
| sigma | A non-negative integer sigma. |
| l | A positive integer l. |
| g | The group element g. |
| x | The group element x = g^d. |
| r | The order r of g. |
| B_ETA | An upper bound on the search space in eta. |
| B_T | An upper bound on the search space in t. |
| verbose | A flag that may be set to True to print intermediary results and status updates when executing the post-processing algorithm. |

## Return value
The logarithm d, or None, if solving for d fails.

