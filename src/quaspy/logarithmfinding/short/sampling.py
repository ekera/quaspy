""" @brief  A module for sampling a frequency pair (j, k) from the distribution
            induced by Ekerå–Håstad's quantum algorithm for finding a given
            short discrete logarithm d in a group of unknown order r. """

import gmpy2;

from gmpy2 import mpz;
from gmpy2 import mpq;
from gmpy2 import mpfr;

from gmpy2 import const_pi as mpfr_const_pi;
from gmpy2 import cos as mpfr_cos;

from ...math.modular import truncmod;

from ...math.random import sample_integer;

from ..sampling import B_DEFAULT_DELTA;
from ..sampling import B_DEFAULT_ETA;

B_DEFAULT = 1000;

def sample_j_k_given_d_r_tau(
  d,
  r,
  m,
  l,
  tau,
  verbose = False,
  extended_result = False):

  """ @brief  Samples a frequency pair (j, k) from the distribution induced by
              Ekerå–Håstad's quantum algorithm for finding a given short
              discrete logarithm d in a group of unknown order r.

      To sample the distribution, j is first picked uniformly at random from
      [0, 2^(m + l)). A pivot is then selected uniformly at random from [0, 1),
      and the optimal frequency k0(j) for k computed.

      For offsets 0, ±1, ±2, .., ±(B - 1) from k0(j) the probabilities
      2^(m+l) * P(j, k = (k0(j) + offset) mod 2^l) of observing (j, k) are then
      subtracted from the pivot, and (j, k) returned as soon as pivot <= 0
      provided that (j, k) is tau-good by Def. 1 in [E23p]. Otherwise, None
      is returned. The bound B is setup as a function of tau so that all
      tau-good pairs (j, k) are included in the search.

      Note that it follows from the analysis in [E20] and [E23p] that j is
      distributed uniformly at random when r >= 2^(m + l) + (2^l - 1) * d.
      This function checks that this requirement is respected.

      [E20] Ekerå, M.: "On post-processing in the quantum algorithm for
                       computing short discrete logarithms".
                       Des., Codes and Cryptogr. 88, pp. 2313–2335 (2020).

      [E23p] Ekerå, M.: "On the success probability of the quantum algorithm for
                         the short DLP". ArXiv 2309.01754 (2023).

      @param d  The discrete logarithm d in [1, r).

      @param r  The order r. Used to check that r >= 2^(m + l) + (2^l - 1) * d.
                May be set to None in which case the check is not performed.

      @param m  A positive integer m such that d < 2^m.

      @param l  An non-negative integer l such that m + l is the length of the
                first control register in the quantum algorithm.

                It is furthermore required that r >= 2^(m + l) + (2^l - 1) * d,
                as this function is based on the analysis in [E20] that imposes
                this requirement so as to simplify the analysis.

      @param tau  The parameter tau. An integer on (0, l].

      @param verbose  A flag that may be set to True to print intermediary
                      results when sampling.

      @param extended_result  A flag that may be set to True to not only return
                              the frequency pair [j, k], but
                              [[j, k], [k0(j), offset]].

      @return   The frequency pair [j, k] sampled if the extended_result flag is
                set to False, or [[j, k], [k0(j), offset]] if the
                extended_result flag is set to True, or None if sampling failed
                because the upper bound B on the offset from the optimal
                frequency k0(j) was reached, or because (j, k) is not tau-good.
  """

  result = sample_j_k_given_d_r(
             d,
             r,
             m,
             l,
             B = (2 ** tau) + 2,
             verbose = verbose,
             extended_result = True);

  if None == result:
    return None;

  [[j, k], _] = result;

  alpha = truncmod(mpz(d * j + (2 ** m) * k), mpz(2 ** (m + l)));
  if not (abs(alpha) < mpz(2 ** (m + tau))):
    return None;

  if extended_result:
    return result;
  else:
    return [j, k];

def sample_j_k_given_d_r(
  d,
  r,
  m,
  l,
  B = B_DEFAULT,
  verbose = False,
  extended_result = False):

  """ @brief  Samples a frequency pair (j, k) from the distribution induced by
              Ekerå–Håstad's quantum algorithm for finding a given short
              discrete logarithm d in a group of unknown order r.

      To sample the distribution, j is first picked uniformly at random from
      [0, 2^(m + l)). A pivot is then selected uniformly at random from [0, 1),
      and the optimal frequency k0(j) for k computed.

      For offsets 0, ±1, ±2, .., ±(B - 1) from k0(j) the probabilities
      2^(m+l) * P(j, k = (k0(j) + offset) mod 2^l) of observing (j, k) are then
      subtracted from the pivot, and (j, k) returned as soon as pivot <= 0.

      Note that it follows from the analysis in [E20] and [E23p] that j is
      distributed uniformly at random when r >= 2^(m + l) + (2^l - 1) * d. This
      function checks that this requirement is respected.

      [E20] Ekerå, M.: "On post-processing in the quantum algorithm for
                       computing short discrete logarithms".
                       Des., Codes and Cryptogr. 88, pp. 2313–2335 (2020).

      [E23p] Ekerå, M.: "On the success probability of the quantum algorithm for
                         the short DLP". ArXiv 2309.01754 (2023).

      @param d  The discrete logarithm d in [1, r).

      @param r  The order r. Used to check that r >= 2^(m + l) + (2^l - 1) * d.
                May be set to None in which case the check is not performed.

      @param m  A positive integer m such that d < 2^m.

      @param l  An non-negative integer l such that m + l is the length of the
                first control register in the quantum algorithm.

                It is furthermore required that r >= 2^(m + l) + (2^l - 1) * d,
                as this function is based on the analysis in [E20] that imposes
                this requirement so as to simplify the analysis.

      @param B  A parameter B that upper-bounds the offset from the optimal
                frequency k0(j) as explained above.

      @param verbose  A flag that may be set to True to print intermediary
                      results when sampling.

      @param extended_result  A flag that may be set to True to not only return
                              the frequency pair [j, k], but
                              [[j, k], [k0(j), offset]].

      @return   The frequency pair [j, k] sampled if the extended_result flag is
                set to False, or [[j, k], [k0(j), offset]] if the
                extended_result flag is set to True, or None if sampling failed
                because the upper bound B on the offset from the optimal
                frequency k0(j) was reached. """

  # Sanity checks.
  if (m <= 0) or (not (d < (2 ** m))):
    raise Exception("Error: Incorrect parameter m.");

  if l <= 0:
    raise Exception("Error: Incorrect parameter l.");

  # Check the requirement that r >= 2 ** (m + l) + (2 ** l - 1) * d.
  if (r != None) and (r < 2 ** (m + l) + (2 ** l - 1) * d):
    raise Exception("Error: It is required that r >= 2^(m + l) + (2^l - 1) * d "
      + "as this function is based on the analysis in [E20] that imposes this "
      + "requirement so as to simplify the analysis.");

  # Sample the frequency j uniformly at random from [0, 2^(m + l)).
  j = mpz(sample_integer(2 ** (m + l)));
  if verbose:
    print("Sampled j =", str(j));

  # Compute the optimal frequency k0(j).
  k0 = mpz(-round(mpq(d * j, 2 ** m))) % mpz(2 ** l);
  if verbose:
    print("Computed k0(j) =", k0);

  # Define precision.
  precision = 2 * (m + 2 * l);
  if precision < 256:
    precision = 256;

  # Swap out the current precision and set the new precision.
  swapped_out_precision = gmpy2.get_context().precision;
  gmpy2.get_context().precision = precision;

  # Explore a region around k0(j).
  pivot = \
    mpfr(mpz(sample_integer(mpz(2 ** precision))), precision) / \
      mpfr(mpz(2 ** precision), precision);

  if verbose:
    print("Sampled pivot =", str(pivot) + "\n");

  def P(j, k):
    # Compute the probability as described in Sect. 3.2 of [E20].
    alpha = truncmod(mpz(d * j + (2 ** m) * k), mpz(2 ** (m + l)));

    if alpha == 0:
      # Use the expression for P(0), see Sect. 3.2 in [E20].
      result  = ((2 ** (m + l)) - ((2 ** l) - 1) * d) * (2 ** (2 * l));
      result += (((2 ** l) * d) / 3) * ((2 ** l) - 1) * ((2 ** (l + 1)) - 1);
      result /= 2 ** (2 * (m + 2 * l));

    else:
      # Use the expression for P(theta), see Sect. 3.2 in [E20].
      theta = mpfr(2 * mpfr_const_pi(precision) * alpha) / mpz(2 ** (m + l));

      result  = mpfr_cos(mpz((2 ** l) - 1) * theta) - \
          mpfr_cos(mpz(2 ** l) * theta);
      result /= 1 - mpfr_cos(theta);
      result -= 1;
      result /= 2;
      result  = mpz((2 ** l) - 1) - result;
      result *= mpz(2 * d);
      result += (mpz(2 ** (m + l)) - mpz((2 ** l) - 1) * d) * \
        (1 - mpfr_cos(mpz(2 ** l) * theta));
      result /= 1 - mpfr_cos(theta);
      result /= mpz(2 ** (2 * (m + 2 * l)));

    return result;

  # Sample k.
  if verbose:
    print("Sampling k...");

  for offset in range(B):
    for sign in [1, -1]:
      if (0 == offset) and (-1 == sign):
        continue;

      k = mpz(k0 + sign * offset) % mpz(2 ** l);
      probability = mpz(2 ** (m + l)) * P(j, k);
      pivot -= probability;

      if verbose:
        print("Offset:", sign * offset, "-- Probability:", probability, \
          "-- Pivot:", pivot, "-- k:", k);

      if pivot <= 0:
        # Restore precision.
        gmpy2.get_context().precision = swapped_out_precision;

        if extended_result:
          return [[j, k], [k0, sign * offset]];
        else:
          return [j, k];

  # Restore precision.
  gmpy2.get_context().precision = swapped_out_precision;

  return None;