""" @brief  A module for sampling a frequency pair (j, k) heuristically from the
            distribution induced by the quantum part of Shor's algorithm for
            finding a discrete logarithm d in a group of known order r, or from
            the distribution induced by the quantum part of Ekerå's variation
            thereof, depending on how parameters are selected. """

from ..sampling import sample_j_k_given_d_r_heuristic as \
  sample_j_k_given_d_r_heuristic_general;

from ...utils.timeout import Timeout;

from ..sampling import B_DEFAULT_DELTA;
from ..sampling import B_DEFAULT_ETA;
from ..sampling import DEFAULT_INTEGRATION_STEPS;

from gmpy2 import mpz;

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
              distribution induced by the quantum part of Shor's algorithm for
              finding a discrete logarithm d in a group of known order r, or
              from the distribution induced by the quantum part of Ekerå's
              variation thereof, depending on how parameters are selected.

      The sampling procedure is described in Sect. 5 of [E19p].

      [E19p] Ekerå, M.: "Revisiting Shor's quantum algorithm for computing
                         general discrete logarithms".
                        ArXiv 1905.09084v4 (2024).

      @param d  The discrete logarithm d in [1, r).

      @param r  The order r.

      @param m  A positive integer m such that r < 2^m.

      @param sigma  A non-negative integer sigma such that m + sigma is the
                    length of the first control register in the quantum part of
                    the algorithm.

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

  return sample_j_k_given_d_r_heuristic_general(
           d = d,
           r = r,
           m = m,
           sigma = sigma,
           l = l,
           B_DELTA = B_DELTA,
           B_DEFAULT_ETA = B_ETA,
           integration_steps = integration_steps,
           verbose = verbose,
           timeout = timeout,
           extended_result = extended_result);
