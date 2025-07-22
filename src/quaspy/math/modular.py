""" @brief  A module for modular arithmetic. """

from gmpy2 import mpz;

def truncmod(x : int | mpz, N : int | mpz) -> int | mpz:

  """ @brief  Returns x mod N constrained to the interval [N/2, N/2).

      @param x  The integer x.

      @param N  The modulus N.

      @return The integer x mod N, constrained to the interval [N/2, N/2). """

  x = x % N;

  if 2 * x >= N:
    x -= N;

  return x;