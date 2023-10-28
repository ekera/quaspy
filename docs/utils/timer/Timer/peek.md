## Method: <code>Timer.peek(self)</code>
Peeks at the timer, returning the number of seconds elapsed.

## Import directive
```python
from quaspy.utils.timer import Timer
```

## Parent class
- [<code>Timer</code>](../Timer.md)

## Prototype
```python
def peek(self)
```

## Return value
If the timer is stopped, the time delta is returned. Otherwise, the sum of the time delta and the current offset of the timer returned.

