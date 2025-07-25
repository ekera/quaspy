## Function: <code>expected\_u\_for\_j\_k\_d(j, k, m, l, d, ..)</code>
Computes the vector u that is associated with a given frequency pair (j, k) and logarithm d, when using the post-processing algorithm from [[E23p]](https://doi.org/10.48550/arXiv.2309.01754) to post-process the output of the quantum part of Ekerå–Håstad's algorithm [[EH17]](https://doi.org/10.1007/978-3-319-59879-6_20).

Recall that the vector v = (truncmod(-2^m k, 2^(m+l)), 0), and that the vector u, which is in the lattice L^tau, is such that the difference u - v = (truncmod(dj - 2^m k, 2^(m+l)), 2^tau d).

Given j, k, m, l, d and tau this function returns the vector u.

## Import directive
```python
from quaspy.logarithmfinding.short.postprocessing import expected_u_for_j_k_d
```

## Parent module
- [<code>postprocessing</code>](README.md)

## Prototype
```python
def expected_u_for_j_k_d(j : int | gmpy2.mpz,
                         k : int | gmpy2.mpz,
                         m : int,
                         l : int,
                         d : int | gmpy2.mpz,
                         tau : int = 0)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| j | The frequency j. An integer on [0, 2^(m + l)). |
| k | The frequency k. An integer on [0, 2^l). |
| m | A positive integer m such that d < 2^m. |
| l | A positive integer l. The control registers in the quantum part of the algorithm are of lengths m + l and l qubits, respectively. |
| d | The discrete logarithm d. |
| tau | An integer tau on [0, l]. Used to scale the basis for the lattice L^tau that is used in the post-processing, and that is generated by the vectors (j, 2^tau) and (2^(m + l), 0). |

## Return value
The vector u.

