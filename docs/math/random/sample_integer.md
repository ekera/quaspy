## Function: <code>sample\_integer(B)</code>
Returns an integer selected uniformly at random from [0, B).

> This function calls randbelow() to select an integer uniformly at random from [0, B). In practice, randbelow() returns an integer that may be conjectured to be indistinguishable from an integer that is selected uniformly at random from [0, B).

## Import directive
```python
from quaspy.math.random import sample_integer
```

## Parent module
- [<code>random</code>](README.md)

## Prototype
```python
def sample_integer(B : int)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| B | The upper bound B. |

## Return value
An integer selected uniformly at random from [0, B).

