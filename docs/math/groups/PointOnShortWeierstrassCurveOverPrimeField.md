## Class: <code>PointOnShortWeierstrassCurveOverPrimeField</code>
A class that represents a point (x, y) on a short Weierstrass curve over a prime field.

Note that we write elliptic curve groups multiplicatively: This is non-standard, but beneficial since it allows elliptic curve groups to be substituted for multiplicative groups, and vice versa, in function calls. More specifically, P + Q, [e] P and [-1] P becomes P * Q, P^e and P^-1, respectively, when writing the group multiplicatively.

## Import directive
```python
from quaspy.math.groups import PointOnShortWeierstrassCurveOverPrimeField
```

## Parent module
- [<code>groups</code>](README.md)

## Methods
- [<code>\_\_eq\_\_(self, Q)</code>](PointOnShortWeierstrassCurveOverPrimeField/__eq__.md)

  Compares this point to another point Q.

- [<code>\_\_hash\_\_(self)</code>](PointOnShortWeierstrassCurveOverPrimeField/__hash__.md)

  Returns the hash digest of this pointgroup element.

- [<code>\_\_init\_\_(self, x, y, E)</code>](PointOnShortWeierstrassCurveOverPrimeField/__init__.md)

  Constructs the point (x, y) on the elliptic curve E.

- [<code>\_\_mul\_\_(self, Q)</code>](PointOnShortWeierstrassCurveOverPrimeField/__mul__.md)

  Returns the point P * Q, for P this point.

- [<code>\_\_pow\_\_(self, e)</code>](PointOnShortWeierstrassCurveOverPrimeField/__pow__.md)

  Returns the point P^e, for P this point.

- [<code>\_\_repr\_\_(self)</code>](PointOnShortWeierstrassCurveOverPrimeField/__repr__.md)

  Returns a string representation of the group element.

- [<code>\_\_str\_\_(self)</code>](PointOnShortWeierstrassCurveOverPrimeField/__str__.md)

  Returns a string representation of the group element.

- [<code>invert(self)</code>](PointOnShortWeierstrassCurveOverPrimeField/invert.md)

  Returns the point P^-1, for P this point.

- [<code>is_identity(self)</code>](PointOnShortWeierstrassCurveOverPrimeField/is_identity.md)

  Returns True if this point is the identity, False otherwise.

