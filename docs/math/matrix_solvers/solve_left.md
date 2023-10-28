## Function: <code>solve\_left(B, o)</code>
Given a 2 x 2 integer matrix B, and an integer row vector o such that c B = o, this function returns the integer row vector c, or None if the equation c B = o has no integer solution

## Import directive
```python
from quaspy.math.matrix_solvers import solve_left
```

## Parent module
- [<code>matrix_solvers</code>](README.md)

## Prototype
```python
def solve_left(B,
               o)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| B | The 2 x 2 integer matrix B = [b1, b2], where b1 = [b_11, b_22] and b2 = [b_21, b_22] are the row vectors that make up B. |
| o | The integer row vector o = [o1, o2]. |

## Return value
The integer row vector c = [c1, c2] such that c B = o, or None if the equation c B = o has no integer solution.

