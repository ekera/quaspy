## Method: <code>Timer.\_\_add\_\_(a, b)</code>
Adds the time deltas of two stopped timers, returning a new timer with said time delta.

The new timer is left in the stopped state until manually started.

If either timer is running, an exception is raised.

## Import directive
```python
from quaspy.utils.timer import Timer
```

## Parent class
- [<code>Timer</code>](../Timer.md)

## Prototype
```python
def __add__(a,
            b)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| a | The first timer. |
| b | The second timer. |

## Return value
A new timer, in the stopped state, with a time delta equal to the sum of the time deltas of the first and second timers.

