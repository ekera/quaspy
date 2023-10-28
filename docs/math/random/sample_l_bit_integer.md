## Function: <code>sample\_l\_bit\_integer(l)</code>
Returns an l-bit integer selected uniformly at random from the set of all such integers.

> This function calls sample_integer() to select an integer uniformly at random from [0, B) for B = 2^(l-1). In practice, sample_integer() returns an integer that may be conjectured to be indistinguishable from an integer that is selected uniformly at random from [0, B).

## Import directive
```python
from quaspy.math.random import sample_l_bit_integer
```

## Parent module
- [<code>random</code>](README.md)

## Prototype
```python
def sample_l_bit_integer(l)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| l | The bit length l of the prime to sample. |

## Return value
An l-bit integer selected uniformly at random from the set of all such integers.

