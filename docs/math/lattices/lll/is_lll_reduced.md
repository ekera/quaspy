## Function: <code>is\_lll\_reduced(B, ..)</code>
Checks if an n x d basis matrix B is delta-LLL-reduced, where LLL is short for Lenstra–Lenstra–Lovász [[LLL82]](https://doi.org/10.1007/BF01457454).

The basis matrix B is represented as a list of lists [b_1, ..., b_n], where b_i = [b_i1, ..., b_id] for i = 1, ..., n are the n rows of A. It is required that B has integer entries of type int or mpz.

## Import directive
```python
from quaspy.math.lattices.lll import is_lll_reduced
```

## Parent module
- [<code>lll</code>](README.md)

## Prototype
```python
def is_lll_reduced(B : list,
                   delta : float = 0.99,
                   gs : list[list[list[int | gmpy2.mpz | gmpy2.mpq | gmpy2.mpfr]]] | None = None,
                   precision : int | None = None)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| B | The basis matrix B. |
| delta | The delta parameter in the Lovász condition. Must be on the interval (1/4, 1]. |
| gs | The list [Bs, M], where Bs is the Gram–Schmidt orthogonalized basis for B, and M is the associated matrix of Gram–Schmidt projection factors, as returned by calling<br><br>[Bs, M] = gram_schmidt(B, precision = precision),<br><br>or None, in which case this function will make the above call. |
| precision | The precision to use when computing the Gram–Schmidt projection factors in M. May be set to None, in which case projection factors are represented as exact quotients.<br><br>Note that this parameter only has an effect if gs is set to None as the Gram–Schmidt orthogonalized basis Bs and the associated matrix M of Gram–Schmidt projection factors are otherwise pre-computed. |

## Return value
True if the basis B is delta-LLL-reduced, False otherwise.

