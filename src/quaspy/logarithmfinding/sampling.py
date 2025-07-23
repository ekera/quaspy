""" @brief  A module for heuristically sampling frequency pairs from the
            distribution induced by the quantum part of Shor's algorithm
            [Shor94] for finding discrete logarithms in groups of known
            order, modified as in [E19p].

    This module may also be used to sample frequency pairs from the distribution
    induced by the quantum parts of Ekerå–Håstad's [EH17] and Ekerå's [E21]
    variations of said algorithm, depending on how parameters are selected.

    [Shor94] Shor, P.W.: "Algorithms for Quantum Computation: Discrete
                          Logarithms and Factoring".
                         In: Proceedings from FOCS '94, pp. 124–134 (1994).

    [EH17] Ekerå, M. and Håstad, J.: "Quantum Algorithms for Computing Short
                                      Discrete Logarithms and Factoring RSA
                                      Integers.". In: PQCrypto 2017.
                                     Springer LNCS 10346, pp. 347–363 (2017).

    [E19p] Ekerå, M.: "Revisiting Shor's quantum algorithm for computing general
                       discrete logarithms".
                      ArXiv 1905.09084v4 (2024).

    [E21] Ekerå, M.: "Quantum algorithms for computing general discrete
                      logarithms and orders with tradeoffs".
                     J. Math. Cryptol. 15(1), pp. 359–407 (2021). """

import gmpy2;

from gmpy2 import mpz;
from gmpy2 import mpfr;

from gmpy2 import const_pi as mpfr_const_pi;
from gmpy2 import sin as mpfr_sin;
from gmpy2 import log as mpfr_log;

from gmpy2 import ceil as mpfr_ceil;
from gmpy2 import rint as mpfr_round;

from gmpy2 import invert as mpz_inv;

from ..math.modular import truncmod;

from ..math.random import sample_integer;

from ..math.kappa import kappa;

from ..utils.timeout import Timeout;

B_DEFAULT_DELTA = 1000;

B_DEFAULT_ETA = 10000;

DEFAULT_INTEGRATION_STEPS = 128;

def sample_j_k_given_d_r_heuristic(
  d : int | mpz,
  r : int | mpz,
  m : int,
  sigma : int,
  l : int,
  B_DELTA : int = B_DEFAULT_DELTA,
  B_ETA : int = B_DEFAULT_ETA,
  integration_steps : int = DEFAULT_INTEGRATION_STEPS,
  timeout : int | None | Timeout = None,
  verbose : bool = False,
  extended_result : bool = False) -> list[int | mpz] | None:

  """ @brief  Samples a frequency pair (j, k) heuristically from the
              distribution induced by the quantum part of Shor's algorithm
              [Shor94] for finding a discrete logarithm d in a group of known
              order r, modified as in [E19p].

      The sampling procedure is described in Sect. 5 in [E19p].

      This function may also be used to sample from the distribution induced by
      the quantum parts of Ekerå–Håstad's [EH17] and Ekerå's [E21] algorithms,
      depending on how parameters are selected, see Sect. 7 and App. B in [E19p]
      for further details. Note also that as for sampling from the quantum part
      of Ekerå–Håstad's [EH17] algorithm Quaspy provides better options than
      going via this function.

      [Shor94] Shor, P.W.: "Algorithms for Quantum Computation: Discrete
                            Logarithms and Factoring".
                           In: Proceedings from FOCS '94, pp. 124–134 (1994).

      [EH17] Ekerå, M. and Håstad, J.: "Quantum Algorithms for Computing Short
                                        Discrete Logarithms and Factoring RSA
                                        Integers.". In: PQCrypto 2017.
                                       Springer LNCS 10346, pp. 347–363 (2017).

      [E19p] Ekerå, M.: "Revisiting Shor's quantum algorithm for computing
                         general discrete logarithms".
                        ArXiv 1905.09084v4 (2024).

      [E21] Ekerå, M.: "Quantum algorithms for computing general discrete
                        logarithms and orders with tradeoffs".
                       J. Math. Cryptol. 15(1), pp. 359–407 (2021).

      @param d  The discrete logarithm d in [1, r).

      @param r  The order r.

      @param m  A positive integer m such that r < 2^m when sampling from the
                distribution induced by Shor's algorithm or Ekerå's algorithm,
                and such that d < 2^m when sampling from the distribution
                induced by Ekerå-Håstad's algorithm.

      @param sigma  A non-negative integer sigma such that m + sigma is the
                    length of the first control register in the quantum
                    part of the algorithm.

      @param l  A positive integer l such that l is the length of the second
                control register in the quantum part of the algorithm.

      @param B_DELTA  A parameter that upper-bounds the offset from the optimal
                      frequency k0(j) when sampling k given j and eta.

      @param B_ETA  A parameter that upper-bounds eta when sampling j and eta.

      @param integration_steps  The number of steps to perform when integrating
                                the probability distribution.

      @param timeout  A timeout after which a TimeoutError will be raised and
                      the sampling procedure aborted.

                      The timeout may be represented as an integer specifying
                      the timeout in seconds, or as an instance of the Timeout
                      class. May be set to None, as is the default, in which
                      case no timeout is enforced.

      @param verbose  A flag that may be set to True to print intermediary
                      results when sampling.

      @param extended_result  A flag that may be set to True to not only return
                              the frequency pair (j, k) as a list [j, k], but
                              [[j, k], [k0(j), offset, eta]].

      @return   The frequency pair [j, k] sampled if the extended_result flag is
                set to False, or [[j, k], [k0(j), offset, eta]] if the
                extended_result flag is set to True, or None if sampling failed
                because the upper bound on the offset from the optimal frequency
                k0(j) or on eta were reached. """

  # Sanity checks.
  if m <= 0:
    raise Exception("Error: Incorrect parameter m.");

  if (r < 2) or (r >= (2 ** m)):
    raise Exception("Error: Incorrect parameter r.");

  if (d < 1) or (d >= r):
    raise Exception("Error: Incorrect parameter d.");

  if l <= 0:
    raise Exception("Error: Incorrect parameter l.");

  if sigma < 0:
    raise Exception("Error: Incorrect parameter sigma.");

  if B_DELTA < 0:
    raise Exception("Error: Incorrect parameter B_DELTA.");

  if B_ETA < 0:
    raise Exception("Error: Incorrect parameter B_ETA.");

  # Initial setup.
  timeout = Timeout.parse(timeout);
  timeout.check();

  # Define the precision.
  precision = 3 * (m + sigma);
  if precision < 256:
    precision = 256;

  with gmpy2.context(gmpy2.get_context()) as context:

    # Set the precision.
    context.precision = precision;


    # Computes the probability of observing alpha_r.
    def P_alpha_r(alpha_r, r, m, sigma, eta):

      # Define the precision.
      precision = 3 * (m + sigma);
      if precision < 256:
        precision = 256;

      with gmpy2.context(gmpy2.get_context()) as context:

        # Set the precision.
        context.precision = precision;

        # Check trivial case.
        if (eta == 0) and (alpha_r == 0):
          return mpfr(1) / mpfr(r);

        # Pre-compute theta_r.
        theta_r = 2 * mpfr_const_pi(precision) * alpha_r / mpz(2 ** (m + sigma));

        # Compute the probability.
        tmp = theta_r - 2 * mpfr_const_pi(precision) * eta;

        # Use that 2 sin^2(x / 2) = 1 - cos(x) for improved stability:
        # probability  = 2 * (1 - mpfr_cos(tmp * mpz(2 ** (m + sigma)) / r));
        tmp2 = mpfr_sin(tmp * mpz(2 ** (m + sigma)) / r / 2);
        probability  = 2 * (2 * tmp2 ** 2);
        probability /= (tmp * tmp);
        probability *= r / mpz(2 ** (2 * (m + sigma)));

        # Return the probability.
        return probability;


    def sample_j_eta(steps = integration_steps):

      # Define the precision.
      precision = 3 * (m + sigma);
      if precision < 256:
        precision = 256;

      with gmpy2.context(gmpy2.get_context()) as context:

        # Set the precision.
        context.precision = precision;

        # Sample a pivot.
        pivot = \
          mpfr(mpz(sample_integer(mpz(2 ** precision))), precision) / \
            mpfr(mpz(2 ** precision), precision);

        log2_r = mpz(mpfr_ceil(mpfr_log(r) / mpfr_log(2)));

        for abs_eta in range(0, B_ETA + 1):

          # Check the timeout.
          timeout.check();

          for sign_eta in [1, -1]:

            if (abs_eta == 0) and (sign_eta == -1):
              continue;

            eta = sign_eta * abs_eta;

            for abs_i in range (0, 100):

              # Check the timeout.
              timeout.check();

              for sgn_i in [1, -1]:

                if (abs_i == 0) and (sgn_i == -1):
                  continue;

                i = log2_r + abs_i * sgn_i;

                if i >= (m + sigma):
                  continue;
                elif i <= 0:
                  continue;

                base = mpz(2 ** (i - 1));
                step = mpz(mpfr_round(mpfr(base) / mpfr(steps)));

                for t in range(0, steps):
                  # Positive angles. Use Simpson's method.
                  p1 = P_alpha_r(base + step * (t    ), r, m, sigma, eta);
                  pm = P_alpha_r(base + step * (t + mpfr(1 / 2)), \
                                   r, m, sigma, eta);
                  p2 = P_alpha_r(base + step * (t + 1), r, m, sigma, eta);

                  probability = step * (p1 + 4 * pm + p2) / 6;
                  pivot -= probability;

                  if pivot <= 0:
                    # Select this alpha_r-interval, and eta.
                    return [base + step * (t    ), \
                            base + step * (t + 1), eta];

                  # Negative angles. Use Simpson's method.
                  p1 = P_alpha_r(-(base + step * (t    )), r, m, sigma, eta);
                  pm = P_alpha_r(-(base + step * (t + mpfr(1 / 2))), \
                                   r, m, sigma, eta);
                  p2 = P_alpha_r(-(base + step * (t + 1)), r, m, sigma, eta);

                  probability = step * (p1 + 4 * pm + p2) / 6;
                  pivot -= probability;

                  if pivot <= 0:
                    # Select this alpha_r-interval, and eta.
                    return [-(base + step * (t + 1)), \
                            -(base + step * (t    )), eta];

            if verbose:
              print("  eta =", eta, "---", "Pivot:", pivot);

        # Failed to sample.
        return None;


    if verbose:
      print("Sampling j...");

    result = sample_j_eta();
    if None == result:
      # Failed to sample.
      return None;

    # Extract the interval in alpha_r and eta.
    [alpha_r_0, alpha_r_1, eta] = result;

    # Sample alpha_r from the interval in alpha_r.
    alpha_r  = alpha_r_1 + sample_integer(alpha_r_1 - alpha_r_0);

    # Ensure that the alpha_r sampled is admissible.
    kappa_r = kappa(r);
    alpha_r -= alpha_r % (2 ** kappa_r);

    # Sample j uniformly at random from all values of j that yield alpha_r.
    t_r = sample_integer(2 ** kappa_r);

    r = mpz(r);
    alpha_r = mpz(alpha_r);

    j = mpz(alpha_r / (2 ** (kappa_r))) * mpz_inv(mpz(r / (2 ** (kappa_r))), \
      mpz(2 ** (m + sigma))) + mpz(2 ** (m + sigma - kappa_r)) * t_r;
    j = j % mpz(2 ** (m + sigma));

    if verbose:
      print("\nSampled j =", int(j), "eta =", eta);

    # Compute the optimal frequency k0(j, eta).
    tmp = d * j - (mpfr(d) / mpfr(r)) * (alpha_r - mpz(2 ** (m + sigma)) * eta);
    k0 = mpz(-mpfr_round(tmp / mpz(2 ** (m + sigma - l)))) % mpz (2 ** l);

    if verbose:
      print("Computed k0(j, eta) =", k0);

    # Sample a pivot.
    pivot = \
      mpfr(mpz(sample_integer(mpz(2 ** precision))), precision) / \
        mpfr(mpz(2 ** precision), precision);

    # Sample k.
    if verbose:
      print("Sampling k...");

    for offset in range(B_DELTA + 1):

      # Check the timeout.
      timeout.check();

      for sign in [1, -1]:
        if (0 == offset) and (-1 == sign):
          continue;

        k = mpz(k0 + sign * offset) % mpz(2 ** l);

        phi = tmp + mpz(2 ** (m + sigma - l)) * k;
        phi = truncmod(phi, mpz(2 ** (m + sigma)));
        phi = (2 * mpfr_const_pi(precision) * phi) / mpz(2 ** (m + sigma));

        # Use that 2 sin^2(x / 2) = 1 - cos(x) for improved stability:
        # probability = (mpfr_cos(phi * mpz(2 ** l)) - 1) / (mpfr_cos(phi) - 1);
        probability  = (mpfr_sin(phi * mpz(2 ** l) / 2) ** 2) / \
          (mpfr_sin(phi / 2) ** 2);
        probability /= mpz(2 ** (2 * l));

        pivot -= probability;

        if verbose:
          print("  Offset:", sign * offset, "-- Probability:", probability, \
            "-- Pivot:", pivot, "-- k:", k);

        if pivot <= 0:
          alpha_d = truncmod(d * j + mpz(2 ** (m + sigma - l)) * k, \
                      mpz(2 ** (m + sigma)));
          if verbose:
            print("\nalpha_d =", alpha_d, mpfr_log(abs(alpha_d)) / mpfr_log(2));

          if extended_result:
            return [[int(j), int(k)], [int(k0), sign * offset, eta]];
          else:
            return [int(j), int(k)];

    # Failed to sample.
    return None;
