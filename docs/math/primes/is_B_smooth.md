## Function: <code>is\_B\_smooth(d, B)</code>
Tests if the integer d is B-smooth.

As in [[E22p]](https://doi.org/10.48550/arXiv.2201.07791), d is said to be B-smooth if d = p1^e1 * .. pk^ek, for q1, .., qk pairwise distinct primes, and e1, .., ek positive integer exponents, if it holds that qi^ei <= B for all i in [1, k].

## Import directive
```python
from quaspy.math.primes import is_B_smooth
```

## Parent module
- [<code>primes</code>](README.md)

## Prototype
```python
def is_B_smooth(d,
                B)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| d | The integer d. |
| B | The upper bound B. |

## Return value
True if d is B-smooth, False otherwise.

