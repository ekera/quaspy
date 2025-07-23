""" @brief  A module for finding the largest power of two that divides an
            integer. """

from gmpy2 import mpz;

def kappa(x : int | mpz) -> int:

  """ @brief  Given an integer x, this function returns t such that x = 2^t * o,
              for o odd.

      @param x  The integer x.

      @return   A non-negative integer t such that x = 2^t * o, for o odd. """

  if x == 0:
    return 0;

  t = 0;

  while (x % 2) == 0:
    t += 1;
    x //= 2;

  return t;