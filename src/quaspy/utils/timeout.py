""" @brief  A module for a utility class for handling timeouts. """

from __future__ import annotations;

from .timer import Timer;

class Timeout:

  """ @brief  A utility class for handling timeouts. """

  def __init__(self, timeout : int | None = None):

    """ @brief  Initializes the timeout to a specific time in seconds.

        @param timeout  The timeout in seconds. May be set to None in which case
                        no timeout is enforced. """

    self.timer = Timer().start();
    self.timeout = timeout;

  def check(self) -> None:

    """ @brief  Raises a TimeoutError if the timeout is elapsed. """

    if self.is_elapsed():
      raise TimeoutError("Error: The timeout has elapsed.");

  def is_elapsed(self) -> bool:

    """ @brief  Returns True if the timeout is elapsed, False otherwise.

        @return   True if the timeout is elapsed, False otherwise. """

    if None == self.timeout:
      return False;

    return self.timer.peek() > self.timeout;

  def is_indefinite(self) -> bool:

    """ @brief  Returns True if the timeout in indefinite, False otherwise.

        @return   True if the timeout is indefinite, False otherwise. """

    return self.timeout == None;

  def __str__(self) -> str:

    """ @brief  Returns a string representation of the timeout.

        @return   A string representation of the timeout. """

    if None == self.timeout:
      return "An indefinite timeout.";
    else:
      s = Timer.format(self.timer.peek());

      if not self.is_elapsed():
        return "A timeout of " + Timer.format(self.timeout) + " with " + s + \
          " elapsed.";
      else:
        return "A timeout of " + Timer.format(self.timeout) + " that is " + \
          "elapsed.";

  def __repr__(self) -> str:

    """ @brief  Returns a string representation of the timeout.

        @return   A string representation of the timeout. """

    return str(self);

  @staticmethod
  def parse(timeout : int | None | Timeout) -> Timeout:

    """ @brief  Parses a timeout that may be an integer, None or an instance of
                Timeout and returns a corresponding instance of Timeout.

        This is a static convenience function.

        @param timeout  The timeout to parse.

        @return   An instance of Timeout representing the timeout. """

    if isinstance(timeout, Timeout):
      return timeout;

    if None == timeout:
      return Timeout(None);

    if isinstance(timeout, int):
      if timeout < 0:
        raise ValueError("Error: Incorrect timeout.");

      return Timeout(timeout);

    raise ValueError("Error: Incorrect timeout.");