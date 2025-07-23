## Function: <code>gram\_schmidt(B, ..)</code>
Returns the Gram–Schmidt orthogonalization Bs of an n x d basis matrix B, and the associated matrix M of Gram–Schmidt projection factors which is such that B = M Bs.

The integer matrix B = [b_1, ..., b_n], where b_i = [b_i1, ..., b_id] for i = 1, ..., n are the n row vectors of B. The matrices Bs and M are represented in analogy with how B is represented.

## Import directive
```python
from quaspy.math.lattices.gram_schmidt import gram_schmidt
```

## Parent module
- [<code>gram_schmidt</code>](README.md)

## Prototype
```python
def gram_schmidt(B : list,
                 precision : int | None = None)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| B | The matrix B. |
| precision | The precision to use when computing the Gram–Schmidt projection factors in M. May be set to None, as is the default, in which case the projection factors are represented as exact quotients. |

## Return value
The pair [Bs, M], where Bs is the Gram–Schmidt orthogonalization of B and M is the matrix of Gram–Schmidt projection factors.

