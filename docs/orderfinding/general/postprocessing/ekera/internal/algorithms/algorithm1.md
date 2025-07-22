## Function: <code>algorithm1(g, r_tilde, m, ..)</code>
Recovers a multiple rp of r, assuming r_tilde is such that r = d * r_tilde where d is cm-smooth.

This function implements Alg. 1 from [[E24]](https://doi.org/10.1145/3655026).

As in [[E24]](https://doi.org/10.1145/3655026), d is said to be cm-smooth if d = p1^e1 * .. pk^ek, for q1, ..., qk pairwise distinct primes, and e1, ..., ek positive integer exponents, if it holds that qi^ei <= cm for all i in [1, k].

## Import directive
```python
from quaspy.orderfinding.general.postprocessing.ekera.internal.algorithms import algorithm1
```

## Parent module
- [<code>algorithms</code>](README.md)

## Prototype
```python
def algorithm1(g : CyclicGroupElement,
               r_tilde : int | gmpy2.mpz,
               m : int,
               c : int = 1,
               timeout : int | None | quaspy.utils.timeout.Timeout = None)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| g | The element g of order r. |
| r_tilde | The integer r_tilde. |
| m | A positive integer m such that r < 2^m. |
| c | A parameter c >= 1 that specifies the maximum size of the missing cm-smooth component d in r = d * r_tilde. |
| timeout | A timeout after which a TimeoutError will be raised and the computation aborted.<br><br>The timeout may be represented as an integer specifying the timeout in seconds, or as an instance of the Timeout class. May be set to None, as is the default, in which case no timeout is enforced. |

## Return value
A multiple rp of r, assuming that r_tilde is such that r = d * r_tilde where d is cm-smooth, and None otherwise.

