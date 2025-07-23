## Function: <code>prime\_power\_product(B)</code>
Returns the product of q^e, as q runs over all primes <= B, for e the largest non-negative integer exponent such that q^e <= B.

## Import directive
```python
from quaspy.math.primes import prime_power_product
```

## Parent module
- [<code>primes</code>](README.md)

## Prototype
```python
def prime_power_product(B : int)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| B | The upper bound B. |

## Return value
The product of q^e, as q runs over all primes <= B, for e the largest non-negative integer exponent such that q^e <= B.

