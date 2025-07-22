""" @brief  A module for solving frequency pairs yielded by the quantum part of
            Shor's algorithm for the general discrete logarithm problem
            [Shor94], modified as in [E19p] and [E24t] (see Sect. 5.5), for the
            logarithm given the order.

    [Shor94] Shor, P.W.: "Algorithms for Quantum Computation: Discrete
                          Logarithms and Factoring".
                         In: Proceedings from FOCS '94, pp. 124–134 (1994).

    [E19p] Ekerå, M.: "Revisiting Shor's quantum algorithm for computing general
                       discrete logarithms".
                      ArXiv 1905.09084v4 (2024).

    [E24t] Ekerå, M.: "On factoring integers, and computing discrete logarithms
                       and orders, quantumly". PhD thesis, KTH Royal Institute
                      of Technology (2024). """

from enum import Enum;

from math import gcd;

from gmpy2 import mpz;

from gmpy2 import invert as mpz_inv;

from ....math.groups import CyclicGroupElement;

from ....math.modular import truncmod;

from ....math.lattices.babai import babai;

from ....math.lattices.cvp import solve_cvp;

from ....math.lattices.lll import lll;
from ....math.lattices.lll import LLL_DEFAULT_DELTA;
from ....math.lattices.lll import LLL_DEFAULT_REDUCED_PRECISION;

from ....math.lattices.enumerate import enumerate as enumerate_lattice;

from ....utils.timeout import Timeout;


B_DEFAULT_ETA = 1000;

B_DEFAULT_T = 1000;


class EnumerationOptions(Enum):

  """ @brief  An enumeration of options for solving a list of n frequency pairs
              [[j_1, k_1], ..., [j_n, k_n]] for d given r.

      The meaning of the options are explained in the documentation for the
      solve_multiple_j_k_for_d_given_r() function. """

  TRUE = True;

  FALSE = False;

  CVP = "CVP";

  BOUNDED_BY_TAU = "BOUNDED_BY_TAU"


def solve_j_k_for_d_given_r(
  j : int | mpz,
  k : int | mpz,
  m : int,
  sigma : int,
  l : int,
  g : CyclicGroupElement,
  x : CyclicGroupElement,
  r : int | mpz,
  B_ETA : int = B_DEFAULT_ETA,
  B_T : int = B_DEFAULT_T,
  timeout : int | None | Timeout = None,
  verbose : bool = False) -> int | mpz | None:

  """ @brief  Attempts to compute the general discrete logarithm d = log_g x
              given a frequency pair (j, k) yielded by the quantum part of
              Shor's algorithm [Shor94], modified as in [E19p], and the order r
              of g, by using the post-processing algorithm from [E19p].

      More specifically, the post-processing algorithm is in Sect. 6 of [E19p].

      [Shor94] Shor, P.W.: "Algorithms for Quantum Computation: Discrete
                            Logarithms and Factoring".
                           In: Proceedings from FOCS '94, pp. 124–134 (1994).

      [E19p] Ekerå, M.: "Revisiting Shor's quantum algorithm for computing
                         general discrete logarithms".
                        ArXiv 1905.09084v4 (2024).

      Note that this function does not implement meet-in-the-middle-techniques,
      although it is noted in [E23p] that it is possible to use such techniques
      to speed up the post-processing also for general discrete logarithms.

      [E23p] Ekerå, M.: "On the success probability of the quantum algorithm for
                         the short DLP". ArXiv 2309.01754v2 (2025).

      @param j  The frequency j. An integer on [0, 2^(m + sigma)).

      @param k  The frequency k. An integer on [0, 2^l).

      @param m  A positive integer m such that r < 2^m.

      @param sigma  A non-negative integer sigma such that m + sigma is the
                    length of the first control register in the quantum part of
                    the algorithm.

      @param l  A positive integer l such that l is the length of the second
                control register in the quantum part of the algorithm.

      @param g  The group element g.

      @param x  The group element x = g^d.

      @param r  The order r of g.

      @param B_ETA  An upper bound on the search space in eta.

      @param B_T  An upper bound on the search space in t.

      @param timeout  A timeout after which a TimeoutError will be raised and
                      the computation aborted.

                      The timeout may be represented as an integer specifying
                      the timeout in seconds, or as an instance of the Timeout
                      class. May be set to None, as is the default, in which
                      case no timeout is enforced.


      @param verbose  A flag that may be set to True to print intermediary
                      results and status updates when executing the
                      post-processing algorithm.

      @return   The logarithm d, or None, if solving for d fails. """

  # Initial setup.
  timeout = Timeout.parse(timeout);
  timeout.check();

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

    # Check the timeout.
    timeout.check();

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

        # Check the timeout.
        timeout.check();

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


def solve_multiple_j_k_for_d_given_r(
  j_k_list : list[list[int | mpz]],
  m : int,
  sigma : int,
  l : int,
  g : CyclicGroupElement,
  x : CyclicGroupElement,
  r : int | mpz,
  tau : int = 0,
  delta : float = LLL_DEFAULT_DELTA,
  precision : int | None = None,
  enumerate : bool | EnumerationOptions = False,
  timeout : int | None | Timeout = None,
  verbose : bool = False) -> int | mpz | None:

  """ @brief  Attempts to compute the general discrete logarithm d = log_g x
              given a list of n frequency pairs [[j_1, k_1], ..., [j_n, k_n]]
              yielded by n independent runs of the quantum part of Shor's
              algorithm [Shor94] for computing general discrete logarithms,
              modified as in [E19p], and the order r of g, by using the
              lattice-based post-processing algorithm described in [E19p]
              (see Sect. 6) and [E24t] (see Sect. 5.5).

      [Shor94] Shor, P.W.: "Algorithms for Quantum Computation: Discrete
                            Logarithms and Factoring".
                           In: Proceedings from FOCS '94, pp. 124–134 (1994).

      [E19p] Ekerå, M.: "Revisiting Shor's quantum algorithm for computing
                         general discrete logarithms".
                        ArXiv 1905.09084v4 (2024).

      [E24t] Ekerå, M.: "On factoring integers, and computing discrete
                         logarithms and orders, quantumly".
                        PhD thesis, KTH Royal Institute of Technology (2024).

      Note that this function does not implement meet-in-the-middle-techniques,
      although it is noted in [E23p] that it is possible to use such techniques
      to speed up the post-processing also for general discrete logarithms.

      [E23p] Ekerå, M.: "On the success probability of the quantum algorithm for
                         the short DLP". ArXiv 2309.01754v2 (2025).

      @param j_k_list   The n frequency pairs [[j_1, k_1], ..., [j_n, k_n]]
                        where j_1, ..., j_n are integers on [0, 2^(m + sigma))
                        and k_1, ..., k_n are integers on [0, 2^l).

      @param m  A positive integer m such that r < 2^m.

      @param sigma  A non-negative integer sigma such that m + sigma is the
                    length of the first control register in the quantum part
                    of the algorithm.

      @param l  A positive integer l ≈ m / s for s a tradeoff factor such that
                l is the length of the second control register in the quantum
                part of the algorithm.

      @param g  The group element g.

      @param x  The group element x = g^d.

      @param r  The order r of g.

      @param tau  A positive integer tau. Used to scale the basis for the
                  lattice L^tau that is used in the post-processing.

      @param delta  The parameter delta to use when delta-LLL-reducing the basis
                    for the lattice L^tau used in the post-processing. Must be
                    on the interval (1/4, 1]. A polynomial runtime in the
                    dimension of the lattice is only guaranteed for delta < 1.

      @param precision  The precision to use when computing the Gram–Schmidt
                        projection factors as a part of delta-LLL-reducing the
                        basis for the lattice L^tau used in the post-processing.

                        The precision may be set to None, as is the default, in
                        which case the projection factors are represented as
                        exact quotients.

      @param enumerate  A flag that may be set to True to enumerate vectors in
                        the lattice L^tau (until d is found or the specified
                        timeout has elapsed), or to EnumerationOptions.CVP to
                        consider only a closest vector in the lattice as
                        returned by performing a limited enumeration, or to
                        False to consider only the vector returned by Babai's
                        nearest plane algorithm.

                        May also be set to EnumerationOptions.BOUNDED_BY_TAU in
                        which case all vectors within distance R of the origin
                        of the lattice L^tau are enumerated, where R depends on
                        tau as R = sqrt(n + 1) * 2^(m + sigma - l + tau).

      @param timeout  A timeout after which a TimeoutError will be raised and
                      the computation aborted.

                      The timeout may be represented as an integer specifying
                      the timeout in seconds, or as an instance of the Timeout
                      class. May be set to None, as is the default, in which
                      case no timeout is enforced.

      @param verbose  A flag that may be set to True to print intermediary
                      results and status updates when executing the
                      post-processing algorithm.

      @return   The logarithm d, or None, if solving for d fails. """

  # Initial setup.
  timeout = Timeout.parse(timeout);
  timeout.check();

  if None in j_k_list:
    return None;

  n = len(j_k_list);
  if n < 1:
    return None;

  js, ks = zip(*j_k_list)

  # Setup a basis for the lattice.
  if verbose:
    print("Setting up the basis...");

  row = [];

  for i in range(n):
    alpha_r_i = mpz(truncmod(r * js[i], 2 ** (m + sigma)));
    row.append((mpz(r * js[i]) - alpha_r_i) * mpz(2 ** l));

  row.append(mpz(2 ** (sigma + tau)) * mpz(r));

  A = [row];

  for i in range(n):
    row = (n + 1) * [mpz(0)];
    row[i] = mpz(2 ** (m + sigma + l)) * mpz(r);
    A.append(row);

  # Reduce the basis for the lattice.
  if verbose:
    print("Computing the LLL-reduced basis with reduced precision...");

  B = lll(A,
          delta = delta,
          timeout = timeout,
          precision = LLL_DEFAULT_REDUCED_PRECISION);
            # Pre-compute with reduced precision.

  # Check the timeout.
  timeout.check();

  if verbose:
    print(" ** Refining the LLL-reduced basis with full precision...");

  [B, [Bs, M]] = lll(B,
                     delta = delta,
                     timeout = timeout,
                     gs = True,
                     precision = precision);

  # Check the timeout.
  timeout.check();

  # Setup v.
  if verbose:
    print("Setting up the target vector v...");

  v = [-(mpz(2 ** (m + sigma)) * mpz(k) * mpz(r)) for k in ks] + [mpz(0)];

  if enumerate in [EnumerationOptions.TRUE, True,
                   EnumerationOptions.FALSE, False,
                   EnumerationOptions.CVP]:

    # Find a vector u in the lattice that is close to v.
    if enumerate == EnumerationOptions.CVP:
      if verbose:
        print("Executing an enumeration to solve CVP...");

      # Let u be a closest vector to v in the lattice.
      u = solve_cvp(B, v, timeout = timeout, gs = [Bs, M]);
    else:
      if verbose:
        print("Executing Babai's nearest plane algorithm...");

      # Let u be a close vector to v in the lattice.
      u = babai(B, v, Bs);

    # Check if the correct solution was found.
    if verbose:
      print(" ** Checking:", u);

    candidate_d = u[-1] // mpz(2 ** (sigma + tau)) // mpz(r);

    if -r < candidate_d < 0:
      # Probably d - r was recovered instead of d.
      candidate_d += r;

    if 0 <= candidate_d < r:
      if g ** candidate_d == x:
        return candidate_d;

    if enumerate not in [True, EnumerationOptions.TRUE]:
      return None;

    # Check the timeout.
    timeout.check();

    if verbose:
      print("Enumerating vectors around v...");

    # Set the enumeration radius R = |u - v|^2.
    R2 = sum([(u[i] - v[i]) ** 2 for i in range(n + 1)]);

    # Enumerate with gradually increasing radii.
    while True:

      # Check the timeout.
      timeout.check();

      # Enumerate, whilst checking the timeout.
      for u in enumerate_lattice(B,
                                 R2,
                                 v,
                                 timeout = timeout,
                                 gs = [Bs, M]):

        # Check if the correct solution was found.
        if verbose:
          print(" ** Checking:", u);

        candidate_d = u[-1] // mpz(2 ** (sigma + tau)) // mpz(r);

        if -r < candidate_d < 0:
          # Probably d - r was recovered instead of d.
          candidate_d += r;

        if 0 <= candidate_d < r:
          if g ** candidate_d == x:
            return candidate_d;

      # Double the enumeration radius.
      if verbose:
          print(" ** Doubling the enumeration radius.");

      R2 = 4 * R2;

  elif enumerate == EnumerationOptions.BOUNDED_BY_TAU:

    # Set the enumeration radius R.
    R2 = mpz(n + 1) * mpz(2 ** (2 * (m + sigma + tau))) * mpz(r * r);

    # Enumerate, whilst checking the timeout.
    for u in enumerate_lattice(B,
                               R2,
                               v,
                               timeout = timeout,
                               gs = [Bs, M]):

      # Check if the correct solution was found.
      if verbose:
        print(" ** Checking:", u);

      candidate_d = u[-1] // mpz(2 ** (sigma + tau)) // mpz(r);

      if -r < candidate_d < 0:
        # Probably d - r was recovered instead of d.
        candidate_d += r;

      if 0 <= candidate_d < r:
        if g ** candidate_d == x:
          return candidate_d;

    # No solution was found.
    return None;

  else:
    raise Exception("Error: Incorrect parameters: Unknown enumerate option.");