## Function: <code>lagrange(A, ..)</code>
Returns the Lagrange-reduced basis for a 2 x 2 basis matrix A.

The basis matrix A is represented as as list [u, v], where u = [u_1, u_2] and v = [v_1, v_2] represent the two row vectors that make up the basis.

## Import directive
```python
from quaspy.math.lattices.lagrange import lagrange
```

## Parent module
- [<code>lagrange</code>](README.md)

## Prototype
```python
def lagrange(A : list,
             multiples = None)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| A | The 2 x 2 basis basis matrix A = [u, v], where u = [u_1, u_2] and v = [v_1, v_2] represent the two row vectors of the basis. |
| multiples | Row multiples of the form [[m_uu, m_uv], [m_vu, m_vv]] such that [m_uu * u + m_uv * v, m_vu * u + m_vv * v] is a close to Lagrange-reduced basis for A, or None as is the default if no such side information is available.<br><br>If row multiples are provided, they are used to start the reduction process with a close to reduced basis, typically speeding up the reduction process.<br><br>The matrix representing the row multiple must have full rank if row multiples are provided. For as long as this basic requirement is met, a Lagrange-reduced basis for A will be returned even if the row multiples are way off. |

## Return value
The pair [A', multiples'], where A' = [u', v'] is a Lagrange-reduced basis for A = [u, v], and multiples' is of the form [[m'_uu, m'_uv], [m'_vu, m'_vv]] and A' = [u', v'] = [m'_uu * u + m'uv * v, m'_vu * u + m'_vv * v].

