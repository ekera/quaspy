""" @brief  A module for solving a frequency pair (j, k) yielded by the quantum
            part of Shor's algorithm for the general discrete logarithm problem
            for the logarithm d given the order r. """

from gmpy2 import mpz;

from gmpy2 import invert as mpz_inv;

from math import gcd;

from ....math.groups import CyclicGroupElement;

from ....math.modular import truncmod;

B_DEFAULT_ETA = 1000;

B_DEFAULT_T = 1000;

def solve_j_k_for_d_given_r(
  j,
  k,
  m,
  sigma,
  l,
  g: CyclicGroupElement,
  x: CyclicGroupElement,
  r,
  B_ETA = B_DEFAULT_ETA,
  B_T = B_DEFAULT_T,
  verbose = False):

  """ @brief  Attempts to compute the general discrete logarithm d given a
              frequency pair (j, k) yielded by the quantum part of
              Shor's algorithm as modified in [E19p], and the order r, by using
              the modified post-processing algorithm described in [E19p].

      The modified post-processing algorithm is described in Sect. 6 of [E19p]:

      [E19p] Ekerå, M.: "Revisiting Shor's quantum algorithm for computing
                         general discrete logarithms".
                        ArXiv 1905.09084v3 (2023).

      Note that this function does not implement meet-in-the-middle-techniques,
      although it is noted in [E23p] that it is possible to use such techniques
      to speed up the post-processing also for general discrete logarithms.

      [E23p] Ekerå, M.: "On the success probability of the quantum algorithm for
                         the short DLP". ArXiv 2309.01754 (2023).

      @param j  The frequency j.

      @param k  The frequency k.

      @param m  A positive integer m such that d < 2^m.

      @param sigma  A non-negative integer sigma.

      @param l  A positive integer l.

      @param g  The group element g.

      @param x  The group element x = g^d.

      @param r  The order r of g.

      @param B_ETA  An upper bound on the search space in eta.

      @param B_T  An upper bound on the search space in t.

      @param verbose  A flag that may be set to True to print intermediary
                      results and status updates when executing the
                      post-processing algorithm.

      @return   The logarithm d, or None, if solving for d fails. """

  j = mpz(j);
  k = mpz(k);

  r = mpz(r);

  # Sanity checks.
  if (m <= 0) or (not (r < (2 ** m))):
    raise Exception("Error: Incorrect parameter m.");

  if l <= 0:
    raise Exception("Error: Incorrect parameter l.");

  if sigma < 0:
    raise Exception("Error: Incorrect parameter sigma.");

  # Pre-compute z and w.
  z = (r * j - truncmod(mpz(r * j), mpz(2 ** (m + sigma)))) // \
    mpz(2 ** (m + sigma));
  w = (r * k - truncmod(mpz(r * k), mpz(2 ** l))) // mpz(2 ** l);

  for abs_eta in range(B_ETA):
    for sgn_eta in [1, -1]:
      if (abs_eta == 0) and (sgn_eta == -1):
        continue;
      eta = sgn_eta * abs_eta;

      if gcd(mpz(z + eta), mpz(r)) != 1:
        continue;

      # Precompute (z + eta)^-1 (mod r).
      z_plus_eta_inv = mpz_inv(mpz(z + eta), mpz(r));

      # Recall that: candidate_d = (t - w) * z_plus_eta_inv
      #
      # As t increases or decreases, we step d by w * z_plus_eta_inv or by
      # -w * z_plus_eta_inv, respectively.
      #
      # Let us pre-compute g to the power of these two exponents:
      x_step_plus  = g ** z_plus_eta_inv;
      x_step_minus = x_step_plus ** -1; # = g ** (-z_plus_eta_inv);

      # The starting point x0 is -w * z_plus_eta_inv.
      #
      # Again, let us pre-compute g to the power of this exponent:
      x0 = x_step_minus ** w; # = g ** (-w * z_plus_eta_inv)

      # Since we search from x0 out, setup storage for intermediary values:
      [candidate_x_plus, candidate_x_minus] = [x0, x0];

      for abs_t in range(B_T):
        for sgn_t in [1, -1]:
          if (abs_t == 0) and (sgn_t == -1):
            continue;
          t = sgn_t * abs_t;

          # Compute candidate_d:
          candidate_d = ((t - w) * z_plus_eta_inv) % r;

          # Set candidate_x = g^candidate_d by using the pre-computed values:
          if sgn_t == 1:
            candidate_x = candidate_x_plus;
            candidate_x_plus = candidate_x_plus * x_step_plus;
          else:
            candidate_x_minus = candidate_x_minus * x_step_minus;
            candidate_x = candidate_x_minus;

          if verbose:
            print("Candidate d:", candidate_d, "t:", t, "eta:", eta);

          if candidate_x == x:
            return candidate_d;

  # Failed to solve for d.
  return None;