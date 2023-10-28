## Function: <code>crt(values, moduli)</code>
Given values = [v1, .., vn] and moduli = [N1, .., Nn], this function returns an integer v in [0, N) such that v = vi (mod Ni) for all i in [1, n], where N = N1 * .. * Nn.

> It is required that all Ni >= 2, and pairwise coprime.

## Import directive
```python
from quaspy.math.crt import crt
```

## Parent module
- [<code>crt</code>](README.md)

## Prototype
```python
def crt(values,
        moduli)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| values | The values = [v1, .., vn]. |
| moduli | The moduli = [N1, .., Nn]. |

## Return value
An integer v in [0, N) such that v = vi (mod Ni) for all i in [1, n], where N = N1 * .. * Nn.

