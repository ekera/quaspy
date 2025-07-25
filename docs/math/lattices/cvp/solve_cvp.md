## Function: <code>solve\_cvp(B, t, ..)</code>
Let L be the lattice generated by the rows of the integer basis matrix B, and let t be a target vector in the span of L. Then this function computes and returns a closest vector to t in L.

The basis matrix B is represented as a list of lists [b_1, ..., b_n], where b_i = [b_i1, ..., b_id] for i = 1, ..., n are the n rows of B.

The vector t is analogously represented as a list [t_1, ..., t_d].

It is assumed that the basis B is delta-LLL-reduced.

It is required that B has integer entries of type int or mpz.

Furthermore, it is required that B has full rank and is square since the enumerate() function that this function calls currently imposes these requirements. Said requirements may be relaxed in the future.

## Import directive
```python
from quaspy.math.lattices.cvp import solve_cvp
```

## Parent module
- [<code>cvp</code>](README.md)

## Prototype
```python
def solve_cvp(B : list,
              t : list,
              timeout : int | None | quaspy.utils.timeout.Timeout = None,
              gs : list[list[list[int | gmpy2.mpz | gmpy2.mpq | gmpy2.mpfr]]] | None = None,
              precision : int | None = None)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| B | The matrix B. |
| t | The vector t. |
| timeout | A timeout after which a TimeoutError will be raised and the computation aborted.<br><br>The timeout may be represented as an integer specifying the timeout in seconds, or as an instance of the Timeout class. May be set to None, as is the default, in which case no timeout is enforced. |
| gs | The list [Bs, M], where Bs is the Gram–Schmidt orthogonalized basis for B, and M is the associated matrix of Gram–Schmidt projection factors, as returned by calling<br><br>[Bs, M] = gram_schmidt(B, precision = precision),<br><br>or None, in which case this function will make the above call.<br><br>If you plan to call this functions repeatedly for the same lattice, then time may be saved by not re-computing Bs and M for each call.<br><br>Note that the basis B is typically LLL-reduced by calling lll() before calling this function. The lll() function can return not only B but also Bs and M (since Bs and M are incrementally computed as a part of the LLL reduction process) allowing you to directly pass them along to this function. |
| precision | The precision to use when computing the Gram–Schmidt projection factors. May be set to None, in which case the projection factors are represented as exact quotients.<br><br>Note that this parameter only has an effect if gs is set to None as the Gram–Schmidt orthogonalized basis Bs and the associated matrix M of Gram–Schmidt projection factors are otherwise pre-computed. |

## Return value
A closest vector to t in L, the lattice generated by the rows of the integer basis matrix B.

