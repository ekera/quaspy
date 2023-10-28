""" @brief  A module for sampling a frequency j from the distribution induced by
            Shor's order-finding algorithm for a given order r. """

import gmpy2;

from gmpy2 import mpz;
from gmpy2 import mpfr;

from gmpy2 import const_pi as mpfr_const_pi;
from gmpy2 import rint as mpfr_round;
from gmpy2 import floor as mpfr_floor;
from gmpy2 import cos as mpfr_cos;
from gmpy2 import gcd;

from ...math.modular import truncmod;

from ...math.random import sample_integer;

B_DEFAULT_SAMPLE = 1000;

def optimal_j_for_z_r(z, r, m, l):

  """ @brief  Computes and returns the optimal frequency j0(z) for z in [0, r).

      @param z  The peak index z in [0, r).

      @param r  The order r.

      @param m  A positive integer m such that r < 2^m.

      @param l  An integer l on [0, m), such that m + l is the length of the
                control register in the quantum order-finding algorithm.

      @return   The optimal frequency j0(z) = round(2^(m + l) / r * z). """

  # Sanity checks.
  if (z < 0) or (z >= r):
    raise Exception("Error: Incorrect parameter z.");

  if (m < 0) or (not (r < (2 ** m))):
    raise Exception("Error: Incorrect parameter m.");

  if (l < 0) or (l > m):
    raise Exception("Error: Incorrect parameter l.");

  # Define precision.
  precision = 2 * (m + l);

  # Swap out the current precision and set the new precision.
  swapped_out_precision = gmpy2.get_context().precision;
  gmpy2.get_context().precision = precision;

  # Compute the optimal frequency j0(z).
  z = mpz(z);
  r = mpz(r);

  j0 = mpz(mpfr_round(mpfr(mpz(2 ** (m + l)) * z) / mpfr(r)));

  # Restore precision.
  gmpy2.get_context().precision = swapped_out_precision;

  return j0;


def sample_j_given_r(
  r,
  m,
  l,
  B = B_DEFAULT_SAMPLE,
  verbose = False,
  extended_result = False):

  """ @brief  Samples a frequency j from the distribution induced by Shor's
              order-finding algorithm for a given order r.

      To sample the distribution, one of the r peaks is first picked uniformly
      at random by picking an index z uniformly at random from [0, r). This
      function furthermore selects a pivot uniformly at random from [0, 1).

      For j0(z) the optimal frequency for this peak, it then subtracts the
      probabilities r * P((j0(z) + offset) mod 2^(m + l)) of frequencies that
      are offset by 0, ±1, ±2, .., ±(B - 1) from the optimal offset from the
      pivot, returning j = (j0(z) + offset) mod 2^(m + l) as soon as pivot <= 0.

      @param r  The order r.

      @param m  A positive integer m such that r < 2^m.

      @param l  An integer l on [0, m), such that m + l is the length of the
                control register in the quantum order-finding algorithm.

      @param B  A parameter B that upper-bounds the offset from the optimal
                frequency j0(z) as explained above.

      @param verbose  A flag that may be set to True to print intermediary
                      results when sampling.

      @param extended_result  A flag that may be set to True to not only return
                              the frequency j, but [j, [z, j0(z), offset]].

      @return   The frequency j sampled if the extended_result flag is set to
                False, or [j, [z, j0(z), offset]] if the extended_result flag is
                set to True, or None if sampling failed because the upper bound
                B on the offset from the optimal frequency j0(z) was reached.
      """

  # Sanity checks.
  if (m < 0) or (not (r < (2 ** m))):
    raise Exception("Error: Incorrect parameter m.");

  if (l < 0) or (l > m):
    raise Exception("Error: Incorrect parameter l.");

  # Define precision.
  precision = 2 * (m + l);
  if precision < 256:
    precision = 256;

  # Swap out the current precision and set the new precision.
  swapped_out_precision = gmpy2.get_context().precision;
  gmpy2.get_context().precision = precision;

  # Sample z uniformly at random from [0, r).
  z = mpz(sample_integer(r));
  if verbose:
    print("Sampled z = " + str(z) + "\n");

  r = mpz(r);
  if verbose:
    print("Note: It holds that gcd(z, r) = " + gcd(z, r) + "\n");

  # Compute the optimal frequency j0(z).
  j0 = mpz(mpfr_round(mpfr(mpz(2 ** (m + l)) * z) / mpfr(r)));

  # Explore a region around j0(z).
  pivot = \
    mpfr(mpz(sample_integer(mpz(2 ** precision))), precision) / \
      mpfr(mpz(2 ** precision), precision);

  if verbose:
    print("Sampled pivot =", str(pivot) + "\n");

  def P(j):
    L = mpz(mpfr_floor(mpfr(2 ** (m + l)) / mpfr(r)));

    beta = mpz(2 ** (m + l)) % r;

    alpha_r = truncmod(mpz(r * j), mpz(2 ** (m + l)));

    if alpha_r == 0:
      result = mpfr((L ** 2) * r + (2 * L + 1) * beta);
    else:
      theta_r = mpfr(2 * mpfr_const_pi(precision) * alpha_r) / \
        mpfr(2 ** (m + l));

      result  = beta * (1 - mpfr_cos(theta_r * (L + 1)));
      result += mpfr(r - beta) * (1 - mpfr_cos(theta_r * L));
      result /= 1 - mpfr_cos(theta_r);

    result /= mpz(2 ** (2 * (m + l)));

    return result;

  for offset in range(B):
    for sign in [1, -1]:
      if (0 == offset) and (-1 == sign):
        continue;

      j = (j0 + sign * offset) % (2 ** (m + l));
      probability = r * P(j);
      pivot -= probability;

      if verbose:
        print("Offset:", sign * offset, "-- Probability:", probability, \
          "-- Pivot:", pivot, "-- j:", j);

      if pivot <= 0:
        # Restore precision.
        gmpy2.get_context().precision = swapped_out_precision;

        if extended_result:
          return [j, [z, j0, sign * offset]];
        else:
          return j;

  # Restore precision.
  gmpy2.get_context().precision = swapped_out_precision;

  return None;