## Module: <code>crt</code>
A module for finding the solution to a set of congruence relations with coprime moduli via the Chinese remainder theorem.

## Import directive
```python
import quaspy.math.crt
```

## Parent module
- [<code>math</code>](../README.md)

## Functions
- [<code>crt(values, moduli)</code>](crt.md)

  Given values = [v1, ..., vn] and moduli = [N1, ..., Nn], this function returns an integer v in [0, N) such that v = vi (mod Ni) for all i in [1, n], where N = N1 * ... * Nn.

