## Function: <code>solve\_left(B, t, ..)</code>
Given a full-rank rational n x n matrix B, and an n-dimensional rational row vector t such that c B = t where c is also a rational n-dimensional row vector, this function returns c.

## Import directive
```python
from quaspy.math.matrices import solve_left
```

## Parent module
- [<code>matrices</code>](README.md)

## Prototype
```python
def solve_left(B : list,
               t : list,
               B_inv : list[list[int | gmpy2.mpz | gmpy2.mpq]] | None = None)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| B | The matrix B = [b_1, ..., b_n], for b_i = [b_i1, ..., b_in] for i = 1, ..., n the n rows of B. The entries b_ij for i, j in 1, ..., n must be of type int, mpz or mpq. |
| t | The n-dimensional rational row vector t = [t_1, ..., t_n]. The entries t_i for i = 1, ..., n must be of type int, mpz or mpq. |
| B_inv | The inverse of the matrix B as computed by calling the function inverse(B), or None in which case this function will call said function to compute the inverse. |

## Return value
The n-dimensional row vector c = [c_1, ..., c_n] such that c B = t.

