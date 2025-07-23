## Function: <code>algorithm4(g, S, m, ..)</code>
Returns a subset Sp of S consisting of all r_tilde in S that are such that d * r_tilde is a positive integer multiple of r, where d is cm-smooth.

This function implements Alg. 4 from [[E24]](https://doi.org/10.1145/3655026).

As in [[E24]](https://doi.org/10.1145/3655026), d is said to be cm-smooth if d = p1^e1 * .. pk^ek, for q1, ..., qk pairwise distinct primes, and e1, ..., ek positive integer exponents, if it holds that qi^ei <= cm for all i in [1, k].

## Import directive
```python
from quaspy.orderfinding.general.postprocessing.ekera.internal.algorithms import algorithm4
```

## Parent module
- [<code>algorithms</code>](README.md)

## Prototype
```python
def algorithm4(g : CyclicGroupElement,
               S : set,
               m : int,
               c : int = 1,
               timeout : int | None | quaspy.utils.timeout.Timeout = None)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| g | The element g of order r. |
| S | A set S of candidates for the integer r_tilde. |
| m | A positive integer m such that r < 2^m. |
| c | A parameter c >= 1 that specifies the maximum size of the missing cm-smooth component d in r = d * r_tilde. |
| timeout | A timeout after which a TimeoutError will be raised and the computation aborted.<br><br>The timeout may be represented as an integer specifying the timeout in seconds, or as an instance of the Timeout class. May be set to None, as is the default, in which case no timeout is enforced. |

## Return value
A subset Sp of S consisting of all r_tilde in S that are such that d * r_tilde is a positive integer multiple of r, where d is cm-smooth.

