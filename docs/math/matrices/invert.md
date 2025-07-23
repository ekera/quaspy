## Function: <code>invert(B)</code>
Inverts a full-rank rational n x n-matrix B.

## Import directive
```python
from quaspy.math.matrices import invert
```

## Parent module
- [<code>matrices</code>](README.md)

## Prototype
```python
def invert(B : list)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| B | The matrix B = [b_1, ..., b_n], for b_i = [b_i1, ..., b_in] for i = 1, ..., n the n rows of B. The entries b_ij must be of type int or mpz for integer B, or of type mpq for rational B. |

## Return value
The inverse of the matrix B.

