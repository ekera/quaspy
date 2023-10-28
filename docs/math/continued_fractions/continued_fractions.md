## Function: <code>continued\_fractions(j, m, l, ..)</code>
Expands j / 2^(m + l) in continued fractions, and returns an ordered list of all denominators < 2^((m + l) / 2), unless an upper bound on the denominator is explicitly specified.

## Import directive
```python
from quaspy.math.continued_fractions import continued_fractions
```

## Parent module
- [<code>continued_fractions</code>](README.md)

## Prototype
```python
def continued_fractions(j,
                        m,
                        l,
                        denominator_bound = None)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| j | The frequency j. An integer on [0, 2^(m+l)). |
| m | A positive integer. |
| l | A non-negative integer. |
| denominator_bound | An upper bound on the denominator. If set to None, as is the default, the bound is taken to be 2^((m + l) / 2). |

## Return value
An ordered list of all denominators < 2^((m + l) / 2), unless an upper bound on the denominator is explicitly specified.

