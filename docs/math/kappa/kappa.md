## Function: <code>kappa(x)</code>
Given an integer x, this function returns t such that x = 2^t * o, for o odd.

## Import directive
```python
from quaspy.math.kappa import kappa
```

## Parent module
- [<code>kappa</code>](README.md)

## Prototype
```python
def kappa(x : int | gmpy2.mpz)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| x | The integer x. |

## Return value
A non-negative integer t such that x = 2^t * o, for o odd.

