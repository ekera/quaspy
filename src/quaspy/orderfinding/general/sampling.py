""" @brief  A module for sampling a frequency j from the distribution induced by
            the quantum part of Shor's order-finding algorithm [Shor94] for a
            given order r. This by using the sampling procedure described in
            [E24t] (see in particular Sects. 5.3.5 and 5.4.3).

    This module can also be used to sample from the distribution induced by
    the quantum part of Seifert's variation [Seifert01] of Shor's order-finding
    algorithm for a given order r, depending on how parameters are selected.

    Throughout this module, the algorithms are as described in [E24t], [E24] and
    [E21]. The notation is also inherited from said works.

    [Shor94] Shor, P.W.: "Algorithms for Quantum Computation: Discrete
                          Logarithms and Factoring".
                         In: Proceedings from FOCS '94, pp. 124–134 (1994).

    [Seifert01] Seifert, J.-P.: "Using fewer qubits in Shor's factorization
                                 algorithm via simultaneous Diophantine
                                 approximation". In: CT-RSA 2001.
                                Springer LNCS 2020, pp. 319–227 (2001).

    [E21] Ekerå, M.: "Quantum algorithms for computing general discrete
                      logarithms and orders with tradeoffs".
                     J. Math. Cryptol. 15(1), pp. 359–407 (2021).

    [E24] Ekerå, M.: "On the success probability of quantum order finding".
                     ACM Trans. Quantum Comput. 5(2):11 (2024).

    [E24t] Ekerå, M.: "On factoring integers, and computing discrete logarithms
                       and orders, quantumly".
                      PhD thesis, KTH Royal Institute of Technology (2024). """

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

from ...utils.timeout import Timeout;

B_DEFAULT_SAMPLE = 10 ** 5;

def optimal_j_for_z_r(
  z : int | mpz,
  r : int | mpz,
  m : int,
  l : int) -> int | mpz:

  """ @brief  Computes and returns the optimal frequency j0(z) for z in [0, r).

      @param z  The peak index z in [0, r).

      @param r  The order r.

      @param m  A positive integer m such that r < 2^m.

      @param l  An integer l on [0, m), such that m + l is the length of the
                control register in the quantum part of the algorithm.

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

  with gmpy2.context(gmpy2.get_context()) as context:

    # Set the precision.
    context.precision = precision;

    # Compute the optimal frequency j0(z).
    z = mpz(z);
    r = mpz(r);

    j0 = mpz(mpfr_round(mpfr(mpz(2 ** (m + l)) * z) / mpfr(r)));

    # Return j0(z).
    return j0;


def sample_j_given_r(
  r : int | mpz,
  m : int,
  l : int,
  B : int = B_DEFAULT_SAMPLE,
  timeout : int | None | Timeout = None,
  verbose : bool = False,
  extended_result : bool = False) -> int | list[int, list[int]] | None:

  """ @brief  Samples a frequency j from the distribution induced by the quantum
              part of Shor's order-finding algorithm [Shor94] for a given order
              r. This by using the sampling procedure described in [E24t] (see
              in particular Sects. 5.3.5 and 5.4.3).

      This fuction can also be used to sample from the distribution induced by
      Seifert's variation [Seifert01] of Shor's order-finding algorithm for a
      given order r, depending on how parameters are selected.

      Throughout this function, the algorithms are as described in [E21], [E24]
      and [E24t]. The notation is also inherited from said works.

      To sample the distribution, one of the r peaks is first picked uniformly
      at random by picking an index z uniformly at random from [0, r). This
      function furthermore selects a pivot uniformly at random from [0, 1).

      For j0(z) the optimal frequency for this peak, it then subtracts the
      probabilities r * P((j0(z) + offset) mod 2^(m + l)) of frequencies that
      are offset by 0, ±1, ±2, ..., ±(B - 1) from the optimal offset from the
      pivot, returning j = (j0(z) + offset) mod 2^(m + l) as soon as pivot <= 0.

      [Shor94] Shor, P.W.: "Algorithms for Quantum Computation: Discrete
                            Logarithms and Factoring".
                           In: Proceedings from FOCS '94, pp. 124–134 (1994).

      [Seifert01] Seifert, J.-P.: "Using fewer qubits in Shor's factorization
                                   algorithm via simultaneous Diophantine
                                   approximation". In: CT-RSA 2001.
                                  Springer LNCS 2020, pp. 319–227 (2001).

      [E21] Ekerå, M.: "Quantum algorithms for computing general discrete
                       logarithms and orders with tradeoffs".
                       J. Math. Cryptol. 15(1), pp. 359–407 (2021).

      [E24] Ekerå, M.: "On the success probability of quantum order finding".
                       ACM Trans. Quantum Comput. 5(2):11 (2024).

      [E24t] Ekerå, M.: "On factoring integers, and computing discrete
                         logarithms and orders, quantumly".
                        PhD thesis, KTH Royal Institute of Technology (2024).

      @param r  The order r.

      @param m  A positive integer m such that r < 2^m.

      @param l  An integer l on [0, m), such that m + l is the length of the
                control register in the quantum part of the algorithm.

      @param B  A parameter B that upper-bounds the offset from the optimal
                frequency j0(z) as explained above.

      @param timeout  A timeout after which a TimeoutError will be raised and
                      the sampling procedure aborted.

                      The timeout may be represented as an integer specifying
                      the timeout in seconds, or as an instance of the Timeout
                      class. May be set to None, as is the default, in which
                      case no timeout is enforced.

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

  # Initial setup.
  timeout = Timeout.parse(timeout);
  timeout.check();

  # Define precision.
  precision = 2 * (m + l);
  if precision < 256:
    precision = 256;

  with gmpy2.context(gmpy2.get_context()) as context:

    # Set the precision.
    context.precision = precision;

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

      with gmpy2.context(gmpy2.get_context()) as context:

        # Set the precision.
        context.precision = precision;

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

      # Check the timeout.
      timeout.check();

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
          if extended_result:
            return [int(j), [int(z), int(j0), sign * offset]];
          else:
            return int(j);

    return None;
