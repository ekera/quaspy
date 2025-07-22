""" @brief  A module for vectors. """

from typing import Any;

def norm2(u):

  """ @brief  Returns the square norm of the n-dimensional vector u.

      @param u  The n-dimensional vector u = [u_1, ..., u_n].

      @return   The square norm of u. """

  return sum([ui * ui for ui in u]);