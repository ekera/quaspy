## Method: <code>Timeout.parse(timeout)</code>
Parses a timeout that may be an integer, None or an instance of Timeout and returns a corresponding instance of Timeout.

This is a static convenience function.

## Import directive
```python
from quaspy.utils.timeout import Timeout
```

## Parent class
- [<code>Timeout</code>](../Timeout.md)

## Prototype
```python
def parse(timeout : int | None | Timeout)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| timeout | The timeout to parse. |

## Return value
An instance of Timeout representing the timeout.

