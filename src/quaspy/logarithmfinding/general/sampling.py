""" @brief  A module for sampling a frequency pair (j, k) heuristically from the
            distribution induced by Shor's quantum algorithm for finding a given
            discrete logarithm d in a group of known order r. """

from ..sampling import sample_j_k_given_d_r_heuristic as \
  sample_j_k_given_d_r_heuristic_general;

from ..sampling import B_DEFAULT_DELTA;
from ..sampling import B_DEFAULT_ETA;

def sample_j_k_given_d_r_heuristic(
  d,
  r,
  m,
  l,
  B_DELTA = B_DEFAULT_DELTA,
  B_ETA = B_DEFAULT_ETA,
  verbose = False,
  extended_result = False):

  """ @brief  Samples a frequency pair (j, k) heuristically from the
              distribution induced by Shor's quantum algorithm for finding a
              given short discrete logarithm d in a group of known order r

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

  return sample_j_k_given_d_r_heuristic_general(
           d = d,
           r = r,
           m = m,
           sigma = l,
           l = l,
           B_DELTA = B_DELTA,
           B_DEFAULT_ETA = B_ETA,
           verbose = verbose,
           extended_result = extended_result);