## Function: <code>sample\_l\_bit\_prime(l)</code>
Returns an l-bit prime selected uniformly at random from the set of all such primes.

> This function calls sample_l_bit_integer() to select l-bit integers uniformly at random from the set of all l-bit integers. In practice, sample_l_bit_integer() returns an l-bit integer that may be conjectured to be indistinguishable from an integer that is selected uniformly at random from said set.

> This function assumes that the is_prime() function performs a deterministic primality test. In practice, this test is likely probabilistic, but the probability of incorrectly classifying an composite as a prime is so small that the test may be conjectured to be indistinguishable from a determinstic test.

## Import directive
```python
from quaspy.math.primes import sample_l_bit_prime
```

## Parent module
- [<code>primes</code>](README.md)

## Prototype
```python
def sample_l_bit_prime(l : int)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| l | The bit length l of the prime to sample. |

## Return value
An l-bit prime selected uniformly at random from the set of all such primes.

