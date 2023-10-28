## Function: <code>gram\_schmidt(A)</code>
Returns the Gram–Schmidt orthogonalization of a 2 x 2 matrix A.

The basis matrix A is represented as a list [a1, a2], where a = [a_11, a_22] and a2 = [a_21, a_22] are the two row vectors that form the basis. It is required that A has integer entries.

## Import directive
```python
from quaspy.math.gram_schmidt import gram_schmidt
```

## Parent module
- [<code>gram_schmidt</code>](README.md)

## Prototype
```python
def gram_schmidt(A)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| A | The 2 x 2 integer basis matrix A = [a1, a2], where a1 = [a_11, a_22] and a2 = [a_21, a_22] are the two row vectors that form the basis. |

## Return value
The Gram–Schmidt orthogonalization of A.

