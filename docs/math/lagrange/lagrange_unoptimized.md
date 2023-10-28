## Function: <code>lagrange\_unoptimized(A)</code>
Returns the Lagrange-reduced basis of a 2 x 2 basis matrix A.

> Compared to lagrange() this function is unoptimized.

The basis matrix A is represented as as list [u1, u2], where u1, u2 are row vectors such that u1 = [a_11, a_12] and u2 = [a_21, a_22].

## Import directive
```python
from quaspy.math.lagrange import lagrange_unoptimized
```

## Parent module
- [<code>lagrange</code>](README.md)

## Prototype
```python
def lagrange_unoptimized(A)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| A | The 2 x 2 basis basis matrix A. |

## Return value
The Lagrange-reduced basis of A.

