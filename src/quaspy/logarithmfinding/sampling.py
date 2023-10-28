""" @brief  A module for sampling a frequency pair (j, k) heuristically from the
            distribution induced by Shor's quantum algorithm for finding a given
            discrete logarithm d in a group of known order r, or from the
            distribution induced by Ekerå–Håstad's and Ekerå's quantum
            algorithms, depending on how parameters are selected. """

import gmpy2;

from gmpy2 import mpz;
from gmpy2 import mpfr;

from gmpy2 import const_pi as mpfr_const_pi;
from gmpy2 import sin as mpfr_sin;
from gmpy2 import cos as mpfr_cos;
from gmpy2 import log as mpfr_log;

from gmpy2 import ceil as mpfr_ceil;
from gmpy2 import rint as mpfr_round;

from gmpy2 import invert as mpz_inv;

from ..math.modular import truncmod;

from ..math.random import sample_integer;

from ..math.kappa import kappa;

B_DEFAULT_DELTA = 1000;

B_DEFAULT_ETA = 10000;

def sample_j_k_given_d_r_heuristic(
  d,
  r,
  m,
  sigma,
  l,
  B_DELTA = B_DEFAULT_DELTA,
  B_ETA = B_DEFAULT_ETA,
  verbose = False,
  extended_result = False):

  """ @brief  Samples a frequency pair (j, k) heuristically from the
              distribution induced by Shor's quantum algorithm for finding a
              given short discrete logarithm d in a group of known order r, or
              from the distribution induced by Ekerå–Håstad's and Ekerå's
              quantum algorithms, depending on how parameters are selected.

      The sampling procedure is described in Sect. 5 of [E19p].

      [E19p] Ekerå, M.: "Revisiting Shor's quantum algorithm for computing
                         general discrete logarithms".
                        ArXiv 1905.09084v3 (2023).

      @param d  The discrete logarithm d in [1, r).

      @param r  The order r.

      @param m  A positive integer m such that r < 2^m when sampling from the
                distribution induced by Shor's algorithm or Ekerå's algorithm,
                and such that d < 2^m when sampling from the distribution
                induced by Ekerå-Håstad's algorithm.

      @param sigma  An non-negative integer sigma such that m + sigma is the
                    length of the first control register in the quantum
                    algorithm. It is required that r < 2^(m + sigma).

      @param l  A positive integer l such that l is the length of the second
                control register in the quantum algorithm.

      @param B_DELTA  A parameter that upper-bounds the offset from the optimal
                      frequency k0(j) when sampling k given j and eta.

      @param B_ETA  A parameter that upper-bounds eta when sampling j and eta.

      @param verbose  A flag that may be set to True to print intermediary
                      results when sampling.

      @param extended_result  A flag that may be set to True to not only return
                              the frequency pair [j, k], but
                              [[j, k], [k0(j), offset, eta]].

      @return   The frequency pair [j, k] sampled if the extended_result flag is
                set to False, or [[j, k], [k0(j), offset, eta]] if the
                extended_result flag is set to True, or None if sampling failed
                because the upper bound on the offset from the optimal frequency
                k0(j) or on eta were reached. """

  # Sanity checks.
  if (m <= 0) or (not (d < (2 ** m))):
    raise Exception("Error: Incorrect parameter m.");

  if l <= 0:
    raise Exception("Error: Incorrect parameter l.");

  if sigma < 0:
    raise Exception("Error: Incorrect parameter sigma.");

  if 2 ** (m + sigma) <= r:
    raise Exception("Error: It is required that 2^(m + sigma) > r.");

  # Define the precision.
  precision = 3 * (m + sigma);
  if precision < 256:
    precision = 256;

  # Swap out the current precision and set the new precision.
  swapped_out_precision = gmpy2.get_context().precision;
  gmpy2.get_context().precision = precision;

  # Computes the probability of observing alpha_r.
  def P_alpha_r(alpha_r, r, m, sigma, eta):
    # Define the precision.
    precision = 3 * (m + sigma);
    if precision < 256:
      precision = 256;

    # Check trivial case.
    if (eta == 0) and (alpha_r == 0):
      return mpfr(1) / mpfr(r);

    # Pre-compute theta_r.
    theta_r = 2 * mpfr_const_pi(precision) * alpha_r / mpz(2 ** (m + sigma));

    # Compute the probability.
    tmp = theta_r - 2 * mpfr_const_pi(precision) * eta;

    # Use that 2 sin^2(x / 2) = 1 - cos(x) for improved stability:
    # probability  = 2 * (1 - mpfr_cos(tmp * mpz(2 ** (m + sigma)) / r));
    probability  = 2 * (2 * mpfr_sin(tmp * mpz(2 ** (m + sigma)) / r / 2) ** 2);
    probability /= (tmp * tmp);
    probability *= r / mpz(2 ** (2 * (m + sigma)));

    # Return the probability.
    return probability;

  def sample_j_eta(steps = 128):
    # Define the precision.
    precision = 3 * (m + sigma);
    if precision < 256:
      precision = 256;

    # Sample a pivot.
    pivot = \
      mpfr(mpz(sample_integer(mpz(2 ** precision))), precision) / \
        mpfr(mpz(2 ** precision), precision);

    log2_r = mpz(mpfr_ceil(mpfr_log(r) / mpfr_log(2)));

    for abs_eta in range(0, B_ETA):
      for sign_eta in [1, -1]:

        if (abs_eta == 0) and (sign_eta == -1):
          continue;

        eta = sign_eta * abs_eta;

        for abs_i in range (0, 100):
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
              # Positive angles.
              p1 = P_alpha_r(base + step * (t    ), r, m, sigma, eta);
              p2 = P_alpha_r(base + step * (t + 1), r, m, sigma, eta);

              probability = step * (p1 + p2) / 2;
              pivot -= probability;

              if pivot <= 0:
                # Select this alpha_r-interval, and eta.
                return [base + step * (t    ), \
                        base + step * (t + 1), eta];

              # Negative angles.
              p1 = P_alpha_r(-(base + step * (t    )), r, m, sigma, eta);
              p2 = P_alpha_r(-(base + step * (t + 1)), r, m, sigma, eta);

              probability = step * (p1 + p2) / 2;
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
    # Restore precision.
    gmpy2.get_context().precision = swapped_out_precision;

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

  for offset in range(B_DELTA):
    for sign in [1, -1]:
      if (0 == offset) and (-1 == sign):
        continue;

      k = mpz(k0 + sign * offset) % mpz(2 ** l);

      phi = tmp + mpz(2 ** (m + sigma - l)) * k;
      phi = truncmod(phi, mpz(2 ** (m + sigma)));
      phi = (2 * mpfr_const_pi(precision) * phi) / mpz(2 ** (m + sigma));

      # Use that 2 sin^2(x / 2) = 1 - cos(x) for improved stability:
      # probability  = (mpfr_cos(phi * mpz(2 ** l)) - 1) / (mpfr_cos(phi) - 1);
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

        # Restore precision.
        gmpy2.get_context().precision = swapped_out_precision;

        if extended_result:
          return [[j, k], [k0, sign * offset, eta]];
        else:
          return [j, k];

  # Restore precision.
  gmpy2.get_context().precision = swapped_out_precision;

  # Failed to sample.
  return None;