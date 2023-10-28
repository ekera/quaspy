## Function: <code>split\_N\_given\_g\_r(g, r, N)</code>
Attempts to split N, given the order r of an element g selected uniformly at random from the multiplicative group of the ring of integers modulo N, where g must be explicitly specified. This by using the original algorithm in [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700).

The algorithm succeeds iff r is even and g^(r/2) != -1 (mod N).

## Import directive
```python
from quaspy.factoring.general.postprocessing.shor import split_N_given_g_r
```

## Parent module
- [<code>shor</code>](README.md)

## Prototype
```python
def split_N_given_g_r(g,
                      r,
                      N)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| g | The element g. |
| r | The order r of g mod N. |
| N | The modulus N. As in [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700), it is required that N is odd, and not a perfect prime power. |

## Return value
A set {p, q} of two non-trivial factors of N such that N = pq, or None if the algorithm failed to split N.

