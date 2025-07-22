## Function: <code>is\_valid\_r\_tilde(r_tilde, m)</code>
Checks if r_tilde is an integer on [1, 2^m).

## Import directive
```python
from quaspy.orderfinding.general.postprocessing.ekera.internal.algorithms import is_valid_r_tilde
```

## Parent module
- [<code>algorithms</code>](README.md)

## Prototype
```python
def is_valid_r_tilde(r_tilde : int | gmpy2.mpz,
                     m : int)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| r_tilde | The r_tilde to check. |
| m | The positive integer m. |

## Return value
True if r_tilde is an integer on [1, 2^m), False otherwise.

