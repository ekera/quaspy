## Function: <code>algorithm2(g, r_tilde, m, ..)</code>
Recovers r, assuming r_tilde is such that r = d * r_tilde where d is cm-smooth.

This function implements Alg. 2 from [[E24]](https://doi.org/10.1145/3655026).

As in [[E24]](https://doi.org/10.1145/3655026), d is said to be cm-smooth if d = p1^e1 * .. pk^ek, for q1, .., qk pairwise distinct primes, and e1, .., ek positive integer exponents, if it holds that qi^ei <= cm for all i in [1, k].

## Import directive
```python
from quaspy.orderfinding.general.postprocessing.ekera.internal.algorithms import algorithm2
```

## Parent module
- [<code>algorithms</code>](README.md)

## Prototype
```python
def algorithm2(g,
               r_tilde,
               m,
               c = 1)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| g | The element g of order r. |
| r_tilde | The integer r_tilde. |
| m | A positive integer m such that r < 2^m. |
| c | A parameter c >= 1 that specifies the maximum size of the missing cm-smooth component d in r = d * r_tilde. |

## Return value
The order r, assuming r_tilde is such that r = d * r_tilde where d is cm-smooth. Otherwise, None or some positive integer multiple of r, is returned.

