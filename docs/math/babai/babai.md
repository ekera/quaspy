## Function: <code>babai(A, v, ..)</code>
Uses Babai's nearest plane algorithm [[Babai86]](https://doi.org/10.1007/BF02579403) to find the vector in the lattice generated L by A that is closest to the vector v.

The matrix A is represented as a list [a1, a2], where a = [a_11, a_22] and a2 = [a_21, a_22] are the two row vectors that form the basis. It is required that A is Lagrange reduced, and that A has integer entries.

## Import directive
```python
from quaspy.math.babai import babai
```

## Parent module
- [<code>babai</code>](README.md)

## Prototype
```python
def babai(A,
          v,
          B = None)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| A | The 2 x 2 Lagrange-reduced integer basis matrix A = [a_1, a_2] for the lattice L, where a_1 = [a_11, a_22] and a_2 = [a_21, a_22] are the two row vectors that form the basis. |
| v | The integer row vector v = [v_1, v_2]. |
| B | The Gram–Schmidt orthogonalization B of A, if available, or None in which case the Gram–Schmidt orthogonalization B is computed by this function. (If you plan on calling this function several times for the same lattice L, then you may wish to pre-compute B and pass B along in each call to this function.) |

## Return value
The vector in the lattice L spanned by A that is closest to v.

