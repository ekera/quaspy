""" @brief  A module for continued fractions. """

from gmpy2 import mpz;
from gmpy2 import mpq;
from gmpy2 import mpfr;

def continued_fractions(
    j : int | mpz,
    m : int,
    l : int,
    denominator_bound : bool = None) -> list[int | mpz]:

  """ @brief  Expands j / 2^(m + l) in continued fractions, and returns an
              ordered list of all denominators < 2^((m + l) / 2), unless an
              upper bound on the denominator is explicitly specified.

      @param j  The frequency j. An integer on [0, 2^(m+l)).

      @param m  A positive integer.

      @param l  A non-negative integer.

      @param denominator_bound  An upper bound on the denominator. If set to
                                None, as is the default, the bound is taken to
                                be 2^((m + l) / 2).

      @return   An ordered list of all denominators < 2^((m + l) / 2), unless an
                upper bound on the denominator is explicitly specified. """

  denominators = [];

  if None == denominator_bound:
    denominator_bound = mpfr(2) ** mpfr((m + l) / 2);

  fraction = mpq(mpz(j), mpz(2 ** (m + l)));

  km1 = mpz(0);
  km2 = mpz(1);

  while True:
    integer_part = fraction.numerator // fraction.denominator; # = floor(f)

    denominator = integer_part * km1 + km2;
    if denominator >= denominator_bound:
      break;

    denominators.append(denominator);

    # Update the recursion.
    km2 = km1;
    km1 = denominator;

    # Update the fraction: The next fraction is 1 / (f - floor(f)).
    fraction = fraction - integer_part;
    if fraction == 0:
      # The best approximation has been reached.
      break;

    fraction = mpz(1) / fraction;

  return denominators;