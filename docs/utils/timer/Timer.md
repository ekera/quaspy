## Class: <code>Timer</code>
A utility class for collecting timing statistics.

## Import directive
```python
from quaspy.utils.timer import Timer
```

## Parent module
- [<code>timer</code>](README.md)

## Methods
- [<code>\_\_add\_\_(a, b)</code>](Timer/__add__.md)

  Adds the time deltas of two stopped timers, returning a new timer with said time delta.

- [<code>\_\_init\_\_(self, ..)</code>](Timer/__init__.md)

  Initializes the timer to a specific time delta. The timer is left in the stopped state until it is manually started.

- [<code>\_\_repr\_\_(self)</code>](Timer/__repr__.md)

  Returns a string representation of the timer.

- [<code>\_\_str\_\_(self)</code>](Timer/__str__.md)

  Returns a string representation of the timer.

- [<code>peek(self)</code>](Timer/peek.md)

  Peeks at the timer, returning the number of seconds elapsed.

- [<code>reset(self, ..)</code>](Timer/reset.md)

  Stops the timer and re-initializes it to a specific time delta.

- [<code>restart(self)</code>](Timer/restart.md)

  Resets the timer and then starts it again.

- [<code>start(self)</code>](Timer/start.md)

  Starts the timer if it is currently stopped.

- [<code>stop(self)</code>](Timer/stop.md)

  Stops the timer if it is currently running.

