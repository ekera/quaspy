## Function: <code>expected\_u\_for\_j\_k\_d(j, k, m, l, d, tau)</code>
Computes the vector u that we seek to find for a given frequency pair (j, k) and logarithm d, when using the post-processing algorithms described in [[E20]](https://doi.org/10.1007/s10623-020-00783-2) and [[E23p]](https://doi.org/10.48550/arXiv.2309.01754).

Recall that the vector v = [truncmod(-2^m k, 2^(m+l)), 0], and that the vector u, which is in the lattice L^tau(j), is such that the difference u - v = [truncmod(dj - 2^m k, 2^(m+l)), 2^tau d].

Given j, k, m, l, d and tau this function returns the vector u. In [[E23p]](https://doi.org/10.48550/arXiv.2309.01754) the parameter tau is variable, whereas tau = 0 in [[E20]](https://doi.org/10.1007/s10623-020-00783-2).

## Import directive
```python
from quaspy.logarithmfinding.short.postprocessing import expected_u_for_j_k_d
```

## Parent module
- [<code>postprocessing</code>](README.md)

## Prototype
```python
def expected_u_for_j_k_d(j : int,
                         k : int,
                         m : int,
                         l : int,
                         d : int,
                         tau : int)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| j | The frequency j. An integer on [0, 2^(m + l)). |
| k | The frequency k. An integer on [0, 2^l). |
| m | A positive integer m such that d < 2^m. |
| l | A positive integer l. |
| d | The discrete logarithm d. |
| tau | The parameter tau. An integer on (0, l]. |

## Return value
The vector u.

