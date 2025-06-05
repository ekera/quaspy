## Function: <code>solve\_j\_k\_for\_d(j, k, m, l, g, x, tau, ..)</code>
Attempts to compute the short discrete logarithm d given a frequency pair (j, k) yielded by the quantum part of Ekerå–Håstad's algorithm, by using the post-processing algorithms described in [[E23p]](https://doi.org/10.48550/arXiv.2309.01754).

This function implements the enumeration procedure in Alg. 1 in [[E23p]](https://doi.org/10.48550/arXiv.2309.01754). It is guaranteed to recover d if (j, k) is a tau-good pair, and if j is such that the lattice L^tau(j) is t-balanced, see Defs. 1–3 in [[E23p]](https://doi.org/10.48550/arXiv.2309.01754).

As shown in Thm. 1 in [[E23p]](https://doi.org/10.48550/arXiv.2309.01754), the probability that (j, k) fulfills these conditions for t and tau is at least

max(0, 1 - f(2^tau)) * max(0, 1 - 2^(Delta - 2(t-1) - tau))

for f(B) = 1 / (B - 1) - 1 / (2 (B - 1)^2) - 1 / (6 (B - 1)^3), and for Delta = m - l on [0, m), for m, l parameters to the quantum algorithm.

Furthermore, as shown in Thm. 1 in [[E23p]](https://doi.org/10.48550/arXiv.2309.01754), at most 2^3 c sqrt(N) group operations must be performed during the enumeration provided that a few group elements are first pre-computed, and provided that there is space to store at most 2^3 sqrt(N) / c integers in a lookup table, for c a positive integer constant, and for N = 2^(Delta + tau + 1) + 2^(tau + t + 2) + 2.

## Import directive
```python
from quaspy.logarithmfinding.short.postprocessing import solve_j_k_for_d
```

## Parent module
- [<code>postprocessing</code>](README.md)

## Prototype
```python
def solve_j_k_for_d(j : int,
                    k : int,
                    m : int,
                    l : int,
                    g : CyclicGroupElement,
                    x : CyclicGroupElement,
                    tau : int,
                    t : int = None,
                    c : int = 1,
                    verbose : bool = False)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| j | The frequency j. An integer on [0, 2^(m + l)). |
| k | The frequency k. An integer on [0, 2^l). |
| m | A positive integer m such that d < 2^m. |
| l | An integer l = m - Delta for Delta an integer on [0, m). |
| g | The group element g. |
| x | The group element x = g^d. |
| tau | The parameter tau. An integer on (0, l]. |
| t | The parameter t. An integer on [0, m). May be set to None, in which case t will be implicitly selected so that the lattice L^tau(j) is t-balanced. If t is not set to None, this function will return None if the lattice L^tau(j) is not t-balanced. |
| c | The constant c. A positive integer. |
| verbose | A flag that may be set to True to print intermediary results and status updates when executing the post-processing algorithm. |

## Return value
The logarithm d, or None, if solving for d fails.

