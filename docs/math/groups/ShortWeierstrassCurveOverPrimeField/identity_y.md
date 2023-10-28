## Method: <code>ShortWeierstrassCurveOverPrimeField.identity\_y(self)</code>
Returns the y coordinate of the point used to represent the point at infinity on this curve.

The point at infinity is represented by (0, self.identity_y()), where the y coordinate is 0 if (0, 0) is not in E, and (0, 1) otherwise.

## Import directive
```python
from quaspy.math.groups import ShortWeierstrassCurveOverPrimeField
```

## Parent class
- [<code>ShortWeierstrassCurveOverPrimeField</code>](../ShortWeierstrassCurveOverPrimeField.md)

## Prototype
```python
def identity_y(self)
```

## Return value
The y coordinate of the point used to represent the point at infinity on this curve.

