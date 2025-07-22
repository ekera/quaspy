## Function: <code>is\_lagrange\_reduced(A)</code>
Returns True if the 2 x 2 basis matrix A is Lagrange-reduced, False otherwise.

The basis matrix A is represented as as list [a1, a2], where a1, a2 are row vectors such that a1 = [a_11, a_12] and a2 = [a_21, a_22].

## Import directive
```python
from quaspy.math.lattices.lagrange import is_lagrange_reduced
```

## Parent module
- [<code>lagrange</code>](README.md)

## Prototype
```python
def is_lagrange_reduced(A : list)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| A | The 2 x 2 basis basis matrix A. |

## Return value
True if A is Lagrange-reduced, False otherwise.

