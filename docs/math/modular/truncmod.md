## Function: <code>truncmod(x, N)</code>
Returns x mod N constrained to the interval [N/2, N/2).

## Import directive
```python
from quaspy.math.modular import truncmod
```

## Parent module
- [<code>modular</code>](README.md)

## Prototype
```python
def truncmod(x : int | gmpy2.mpz,
             N : int | gmpy2.mpz)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| x | The integer x. |
| N | The modulus N. |

## Return value
The integer x mod N, constrained to the interval [N/2, N/2).

