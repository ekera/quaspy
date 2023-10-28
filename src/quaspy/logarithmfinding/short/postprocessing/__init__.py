""" @brief  A module for solving a frequency pair (j, k) yielded by the quantum
            part of Ekerå–Håstad's quantum algorithm for the order r. This by
            using the post-processing algorithms in [E20] and [E23p].

    [E20] Ekerå, M.: "On post-processing in the quantum algorithm for computing
                      short discrete logarithms".
                     Des., Codes and Cryptogr. 88, pp. 2313–2335 (2020).

    [E23p] Ekerå, M.: "On the success probability of the quantum algorithm for
                       the short DLP". ArXiv 2309.01754 (2023).

    This implementation does not currently support tradeoffs, as it uses
    Lagrange's lattice basis reduction algorithm. """

from gmpy2 import mpz;
from gmpy2 import mpfr;
from gmpy2 import mpq;

from gmpy2 import floor as mpfr_floor;
from gmpy2 import rint as mpfr_round;
from gmpy2 import ceil as mpfr_ceil;
from gmpy2 import div as mpfr_div;

from gmpy2 import sqrt as mpfr_sqrt;

from ....math.lagrange import lagrange;

from ....math.groups import CyclicGroupElement;

from ....math.modular import truncmod;

from ....math.babai import babai;

from ....math.matrix_solvers import solve_left;

def expected_u_for_j_k_d(
  j : int,
  k : int,
  m : int,
  l : int,
  d: int,
  tau: int):

  """ @brief  Computes the vector u that we seek to find for a given frequency
              pair (j, k) and logarithm d, when using the post-processing
              algorithms described in [E20] and [E23p].

      [E20] Ekerå, M.: "On post-processing in the quantum algorithm for
                        computing short discrete logarithms".
                       Des., Codes and Cryptogr. 88, pp. 2313–2335 (2020).

      [E23p] Ekerå, M.: "On the success probability of the quantum algorithm for
                         the short DLP". ArXiv 2309.01754 (2023).

      Recall that the vector v = [truncmod(-2^m k, 2^(m+l)), 0], and that the
      vector u, which is in the lattice L^tau(j), is such that the difference
      u - v = [truncmod(dj - 2^m k, 2^(m+l)), 2^tau d].

      Given j, k, m, l, d and tau this function returns the vector u. In [E23p]
      the parameter tau is variable, whereas tau = 0 in [E20].

      @param j  The frequency j. An integer on [0, 2^(m + l)).

      @param k  The frequency k. An integer on [0, 2^l).

      @param m  A positive integer m such that d < 2^m.

      @param l  A positive integer l.

      @param d  The discrete logarithm d.

      @param tau  The parameter tau. An integer on (0, l].

      @return   The vector u. """

  # Setup the target vector v.
  v = [truncmod(-(mpz(2 ** m) * k), mpz(2 ** (m + l))), 0];

  # Compute u.
  t1 = mpz(d * j - v[0]);
  t2 = truncmod(t1, mpz(2 ** (m + l)));
  mp = mpz(round(mpq(t1 - t2, mpz(2 ** (m + l)))));

  u = [d * j - mpz(2 ** (m + l)) * mp, mpz(2 ** tau) * mpz(d)];

  # Return u.
  return u;

def solve_j_k_for_d(
  j : int,
  k : int,
  m : int,
  l : int,
  g: CyclicGroupElement,
  x: CyclicGroupElement,
  tau : int,
  t : int = None,
  c : int = 1,
  verbose : bool = False):

  """ @brief  Attempts to compute the short discrete logarithm d given a
              frequency pair (j, k) yielded by the quantum part of
              Ekerå–Håstad's algorithm, by using the post-processing algorithms
              described in [E23p].

      [E23p] Ekerå, M.: "On the success probability of the quantum algorithm for
                         the short DLP". ArXiv 2309.01754 (2023).

      This function implements the enumeration procedure in Alg. 1 in [E23p]. It
      is guaranteed to recover d if (j, k) is a tau-good pair, and if j is such
      that the lattice L^tau(j) is t-balanced, see Def. 1–3 in [E23p].

      As shown in Thm. 1 in [E23p], the probability that (j, k) fulfills these
      conditions for t and tau is at least

          max(0, 1 - f(2^tau)) * max(0, 1 - 2^(Delta - 2(t-1) - tau))

      for f(B) = 1 / (B - 1) - 1 / (2 (B - 1)^2) - 1 / (6 (B - 1)^3), and for
      Delta = m - l on [0, m), for m, l parameters to the quantum algorithm.

      Furthermore, as shown in Thm. 1 in [E23p], at most 2^3 c sqrt(N) group
      operations must be performed during the enumeration provided that a few
      group elements are first pre-computed, and provided that there is space to
      store at most 2^3 sqrt(N) / c integers in a lookup table, for c a positive
      integer constant, and for N = 2^(Delta + tau + 1) + 5 * 2^(tau + t) + 3.

      @param j  The frequency j. An integer on [0, 2^(m + l)).

      @param k  The frequency k. An integer on [0, 2^l).

      @param m  A positive integer m such that d < 2^m.

      @param l  An integer l = m - Delta for Delta an integer on [0, m).

      @param g  The group element g.

      @param x  The group element x = g^d.

      @param tau  The parameter tau. An integer on (0, l].

      @param t  The parameter t. An integer on [0, m). May be set to None, in
                which case t will be implicitly selected so that the lattice
                L^tau(j) is t-balanced. If t is not set to None, this function
                will return None if the lattice L^tau(j) is not t-balanced.

      @param c  The constant c. A positive integer.

      @param verbose  A flag that may be set to True to print intermediary
                      results and status updates when executing the
                      post-processing algorithm.

      @return   The logarithm d, or None, if solving for d fails. """

  # Sanity checks.
  if not (0 < l <= m):
    raise Exception("Error: Incorrect parameter:",
                      "l must be an integer on (0, m].");

  if not (0 <= tau <= l):
    raise Exception("Error: Incorrect parameter:",
                      "tau must be an integer on [0, l].");

  if not (0 <= t < m):
    raise Exception("Error: Incorrect parameter:",
                      "t must be an integer on [0, m).");

  if not (c >= 1):
    raise Exception("Error: Incorrect parameter:",
                      "c must be a positive integer.");

  if not (0 <= j < (2 ** (m + l))):
    raise Exception("Error: Incorrect parameter:",
                      "j must be an integer on [0, 2^(m + l)).");

  if not (0 <= k < (2 ** l)):
    raise Exception("Error: Incorrect parameter:",
                      "k must be an integer on [0, 2^l).");

  j = mpz(j);
  k = mpz(k);

  # Step 1:
  def meet_in_the_middle(g, x, nu1, nu2, B1, B2, s1, s2, c = 1):

    # Step 1.1:
    if B1 < B2:

      # Step 1.1.1:
      return meet_in_the_middle(g, x, nu2, nu1, B2, B1, s2, s1, c);

    # Step 1.2:
    g1 = g ** s1;
    g2 = g ** s2;
    w = (g1 ** nu1) * (g2 ** nu2) * (x ** -1);

    # Step 1.3:
    if B1 == 0:

      if verbose:
        print("Note: Solving without enumerating.")

      # Step 1.3.1:
      if w == g ** 0:

        # Step 1.3.1.1:
        d = nu1 * s1 + nu2 * s2;
        return d;

    # Step 1.4:
    n = c * mpz(mpfr_round(mpfr_sqrt(B1 / (B2 + 1))));

    if verbose:
      print("Computed n =", n);
      print("");

    # Step 1.5:
    if verbose:
      print("The first stage begins:", \
              2 * (ceil(B1 / n) - 1), "operation(s)");

    T = dict();
    T[g ** 0] = 0;

    # Step 1.6:
    s = g1 ** n;
    z_plus = s;
    s_inv = s ** -1;
    z_minus = s_inv;
    i = 1;

    # Step 1.7:
    while True:

      # Step 1.7.1:
      T[z_plus] = i;
      T[z_minus] = -i;

      # Step 1.7.2:
      i = i + 1;
      if i > mpfr_ceil(B1 / n):
        # Step 1.7.2.1:
        break;

      # Step 1.7.3:
      z_plus = z_plus * s;
      z_minus = z_minus * s_inv; # Note: s^-1 was pre-computed above.

    # Step 1.8:
    if verbose:
      print("The second stage begins:", \
              2 * B2 + 2 * (B2 + 1) * (n - 1), "operation(s)");

    z_plus = w;
    z_minus = w;
    j = 0;

    g2_inv = g2 ** -1;

    # Step 1.9:
    while True:

      # Step 1.9.1:
      zp_plus = z_plus;
      zp_minus = z_minus;
      i = 0;

      # Step 1.9.2:
      while True:

        # Step 1.9.2.1:
        if zp_plus in T.keys():
          k = T[zp_plus];

          # Step 1.9.2.1.1:
          d = (nu1 + i - k * n) * s1 + (nu2 + j) * s2;
          return d;

        # Step 1.9.2.2:
        if (j > 0) and (zp_minus in T.keys()):
          k = T[zp_minus];

          # Step 1.9.2.2.1:
          d = (nu1 + i - k * n) * s1 + (nu2 - j) * s2;
          return d;

        # Step 1.9.2.3:
        i = i + 1;
        if i >= n:
          # Step 1.9.2.3.1:
          break;

        # Step 1.9.2.4:
        zp_plus = zp_plus * g1;
        zp_minus = zp_minus * g1;

      # Step 1.9.3:
      j = j + 1;
      if j > B2:
        # Step 1.9.3.1:
        break;

      # Step 1.9.4:
      z_plus = z_plus * g2;
      z_minus = z_minus * g2_inv; # Note: g2^-1 was pre-computed above.

    # Step 1.10:
    return None;

  # Step 2: Setup the basis for the lattice.
  A = [[j, mpz(2 ** tau)], [mpz(2 ** (m + l)), mpz(0)]];

  # Compute the reduced basis for the lattice and extract s1 and s2.
  [B, _] = lagrange(A);

  s1 = B[0];
  s2 = B[1];

  # Compute the norm lambda1 of s1.
  lambda1 = mpfr_sqrt((s1[0] ** 2) + (s1[1] ** 2));

  # Verify the requirement on the norm.
  if None != t:
    if lambda1 < mpz(2 ** (m - t)):
      return None;

  # Computes the inner product between two vectors a and b.
  def inner_product(a, b):
    return a[0] * b[0] + a[1] * b[1];

  # Compute s2_parallel and then s2_perp given s2_parallel.
  mu = mpfr(inner_product(s1, s2)) / mpfr((s1[0] ** 2) + (s1[1] ** 2));

  s2_parallel = [mu * s1[0], mu * s1[1]];
  s2_perp = [s2[0] - s2_parallel[0], s2[1] - s2_parallel[1]];

  # Compute the norm lambda2_perp of s2_perp.
  lambda2_perp = mpfr_sqrt((s2_perp[0] ** 2) + (s2_perp[1] ** 2));

  # Step 4:

  # Define the vector v.
  v = [truncmod(-(2 ** m) * k, 2 ** (m + l)), mpz(0)];

  # Use Babai's algorithm to find the vector o in the lattice closest to v.
  o = babai(B, v);

  # Let nu1 and nu2 be integers such that o = nu1 * s1 + nu2 * s2.
  nu = solve_left(B, o);

  # Sanity check.
  if o != [nu[0] * B[0][0] + nu[1] * B[1][0],
           nu[0] * B[0][1] + nu[1] * B[1][1]]:
    raise Exception("Error: Failed to solve for nu.")

  # Step 5:
  B1 = mpz(mpfr_floor((2 ** (m + tau)) * mpfr_sqrt(mpfr(2)) /
                        mpfr(lambda1) + 1/2));
  B2 = mpz(mpfr_floor((2 ** (m + tau)) * mpfr_sqrt(mpfr(2)) /
                        mpfr(lambda2_perp) + 1/2));

  if verbose:
    print("Computed B1 =", B1);
    print("Computed B2 =", B2);
    print("");

  # Step 6:
  return meet_in_the_middle(g, x,
                            nu[0], nu[1],
                            B1, B2,
                            mpz(s1[1] // (2 ** tau)),
                            mpz(s2[1] // (2 ** tau)),
                            c = c);
