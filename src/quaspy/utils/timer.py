""" @brief  A module for a utility class for collecting timing statistics. """

from __future__ import annotations;

from time import time;
from math import floor;

class Timer:

  """ @brief  A utility class for collecting timing statistics. """

  # A constant for state management indicating that the timer is stopped.
  STOPPED = 0;

  # A constant for state management indicating that the timer is running.
  RUNNING = 1;

  def __init__(self, delta_t = 0):

    """ @brief  Initializes the timer to a specific time delta. The timer is
                left in the stopped state until it is manually started.

        @param delta_t  The time delta. Defaults to zero. """

    self.state = Timer.STOPPED;
    self.delta_t = delta_t;

  def start(self) -> Timer:

    """ @brief  Starts the timer if it is currently stopped.

        @return   The timer. """

    if self.state == Timer.STOPPED:
      self.state = Timer.RUNNING;
      self.t = time();

    return self;

  def stop(self) -> Timer:

    """ @brief  Stops the timer if it is currently running.

        @return   The timer. """

    if self.state == Timer.RUNNING:
      self.state = Timer.STOPPED;
      self.delta_t += time() - self.t;

    return self;

  def reset(self, delta_t = 0) -> Timer:

    """ @brief  Stops the timer and re-initializes it to a specific time delta.

        @param delta_t  The time delta. Defaults to zero.

        @return The timer. """

    self.state = Timer.STOPPED;
    self.delta_t = delta_t;

    return self;

  def restart(self) -> Timer:

    """ @brief  Resets the timer and then starts it again.

        @return   The timer. """

    self.reset();
    self.start();

    return self;

  def peek(self) -> int | float:

    """ @brief  Peeks at the timer, returning the number of seconds elapsed.

        @return   If the timer is stopped, the time delta is returned.
                  Otherwise, the sum of the time delta and the current offset of
                  the timer returned. """

    tmp_delta_t = self.delta_t;
    if self.state == Timer.RUNNING:
      tmp_delta_t += time() - self.t;

    return tmp_delta_t;

  def __add__(a : Timer, b : Timer) -> Timer:

    """ @brief  Adds the time deltas of two stopped timers, returning a new
                timer with said time delta.

        The new timer is left in the stopped state until manually started.

        If either timer is running, an exception is raised.

        @param a   The first timer.

        @param b   The second timer.

        @return   A new timer, in the stopped state, with a time delta equal to
                  the sum of the time deltas of the first and second timers. """

    if (Timer.STOPPED != a.state) or (Timer.STOPPED != b.state):
      raise Exception("Error: Cannot add running timers.");

    return Timer(a.delta_t + b.delta_t);

  def __str__(self) -> str:

    """ @brief  Returns a string representation of the timer.

        @return   A string representation of the timer. """

    return Timer.format(self.peek());

  def __repr__(self) -> str:

    """ @brief  Returns a string representation of the timer.

        @return   A string representation of the timer. """

    return str(self);

  @staticmethod
  def format(delta_t : int | float) -> str:

    """ @brief  Formats a time delta in seconds as a string.

        @param delta_t  The time delta in seconds.

        @return   A string representation of the time delta. """

    # Compute days, hours, minutes, seconds, milliseconds and microseconds.
    days = floor(delta_t / 3600 / 24);
    hours = floor(delta_t / 3600) % 24;
    mins = floor(delta_t / 60) % 60;
    secs = floor(delta_t) % 60;
    ms = int(floor((10 ** 3) * delta_t)) % (10 ** 3);
    us = int(floor((10 ** 6) * delta_t)) % (10 ** 3);

    # Format the time delta as a human-readable string.
    hr = "";

    if days > 0:
      if hr != "":
        hr += " "
      hr += str(days) + " d";

    if hours > 0:
      if hr != "":
        hr += " "
      hr += str(hours) + " h";

    if mins > 0:
      if hr != "":
        hr += " "
      hr += str(mins) + " m";

    if secs > 0:
      if hr != "":
        hr += " "
      hr += str(secs) + " s";

    if ms > 0:
      if hr != "":
        hr += " "
      hr += str(ms) + " ms";

    if us > 0:
      if hr != "":
        hr += " "
      hr += str(us) + " µs";

    # Handle the special case of the time delta being zero.
    if hr == "":
      if isinstance(delta_t, int):
        hr = "0 s";
      else:
        hr = "0 µs";

    return hr;