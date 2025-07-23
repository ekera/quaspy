## Function: <code>lll(A, ..)</code>
Returns the delta-LLL-reduced basis for an n x d basis matrix A, where LLL is short for Lenstra–Lenstra–Lovász [[LLL82]](https://doi.org/10.1007/BF01457454).

The basis matrix A is represented as a list of lists [a_1, ..., a_n], where a_i = [a_i1, ..., a_id] for i = 1, ..., n are the n rows of A. It is required that A has integer entries of type int or mpz.

## Import directive
```python
from quaspy.math.lattices.lll import lll
```

## Parent module
- [<code>lll</code>](README.md)

## Prototype
```python
def lll(A : list,
        delta : float = 0.99,
        timeout : int | None | quaspy.utils.timeout.Timeout = None,
        gs : bool = False,
        precision : int | None = None)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| A | The basis matrix A. |
| delta | The delta parameter in the Lovász condition. Must be on the interval (1/4, 1]. A polynomial runtime in the dimension n of the lattice is only guaranteed for delta < 1. |
| timeout | A timeout after which a TimeoutError will be raised and the computation aborted.<br><br>The timeout may be represented as an integer specifying the timeout in seconds, or as an instance of the Timeout class. May be set to None, as is the default, in which case no timeout is enforced. |
| gs | A flag that may be set to True to return not only the delta-LLL-reduced basis B of A, but also the Gram–Schmidt orthogonalization Bs of B, and the matrix M of Gram–Schmidt projection factors such that B = M Bs. Note that Bs and M are always computed as a part of the LLL reduction process. |
| precision | The precision to use when computing the Gram–Schmidt projection factors in M. May be set to None, as is the default, in which case the projection factors are represented as exact quotients. |

## Return value
The delta-LLL-reduced basis B = [b_1, ..., b_n], where b_i = [b_i1, ..., b_id] for i = 1, ..., n represent the n row vectors that make up the basis, if gs is set to False.

Otherwise, if gs is set to True, [B, [Bs, M]] is returned where B is as above, Bs is the Gram–Schmidt orthogonalization of B and M is the matrix of Gram–Schmidt projection factors. It holds that B = M Bs. The matrices Bs and M are represented in analogy with how A and B are represented.

