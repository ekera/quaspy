## Class: <code>Timeout</code>
A utility class for handling timeouts.

## Import directive
```python
from quaspy.utils.timeout import Timeout
```

## Parent module
- [<code>timeout</code>](README.md)

## Methods
- [<code>\_\_init\_\_(self, ..)</code>](Timeout/__init__.md)

  Initializes the timeout to a specific time in seconds.

- [<code>\_\_repr\_\_(self)</code>](Timeout/__repr__.md)

  Returns a string representation of the timeout.

- [<code>\_\_str\_\_(self)</code>](Timeout/__str__.md)

  Returns a string representation of the timeout.

- [<code>check(self)</code>](Timeout/check.md)

  Raises a TimeoutError if the timeout is elapsed.

- [<code>is_elapsed(self)</code>](Timeout/is_elapsed.md)

  Returns True if the timeout is elapsed, False otherwise.

- [<code>is_indefinite(self)</code>](Timeout/is_indefinite.md)

  Returns True if the timeout is indefinite, False otherwise.

- [<code>parse(timeout)</code>](Timeout/parse.md)

  Parses a timeout that may be an integer, None or an instance of Timeout and returns a corresponding instance of Timeout.

