""" @brief  A module for implementing Algorithms 1–4 from [E22p].

    [E22p] Ekerå, M.: "On the success probability of quantum order finding".
                      ArXiv 2201.07791 (2022). """

from gmpy2 import mpz;

from math import floor;
from math import log;
from math import prod;

from ......math.primes import prime_range;

def is_valid_r_tilde(r_tilde, m):

  """ @brief  Checks if r_tilde is an integer on [1, 2^m).

      @param r_tilde  The r_tilde to check.

      @param m  The positive integer m.

      @return   True if r_tilde is an integer on [1, 2^m), False otherwise. """

  if not (isinstance(r_tilde, int) or isinstance(r_tilde, mpz)):
    return False;

  return (r_tilde >= 1) and (r_tilde < 2 ** m);


def algorithm1(g, r_tilde, m, c = 1):

  """ @brief  Recovers a multiple rp of r, assuming r_tilde is such that
              r = d * r_tilde where d is cm-smooth.

      This function implements Algorithm 1 from [E22p].

      [E22p] Ekerå, M.: "On the success probability of quantum order finding".
                        ArXiv 2201.07791 (2022).

      As in [E22p], d is said to be cm-smooth if d = p1^e1 * .. pk^ek, for
      q1, .., qk pairwise distinct primes, and e1, .., ek positive integer
      exponents, if it holds that qi^ei <= cm for all i in [1, k].

      @param g  The element g of order r.

      @param r_tilde  The integer r_tilde.

      @param m  A positive integer m such that r < 2^m.

      @param c  A parameter c >= 1 that specifies the maximum size of the
                missing cm-smooth component d in r = d * r_tilde.

      @return   A multiple rp of r, assuming that r_tilde is such that
                r = d * r_tilde where d is cm-smooth, and None otherwise. """

  # Step 1:
  if not is_valid_r_tilde(r_tilde, m):
    # Step 1.1:
    return None;

  # Step 2:
  rp = r_tilde; x = g ** r_tilde;

  # Step 3:
  for q in prime_range(floor(c * m) + 1):
    # Step 3.1:
    if x == 1:
      # Step 3.1.1:
      return rp;

    # Step 3.2:
    e = floor(log(c * m) / log(q));

    # Step 3.3:
    q_pow_e = q ** e;
    x = x ** q_pow_e;
    rp = rp * q_pow_e;

  # Step 4:
  if x != 1:
    # Step 4.1:
    return None;

  # Step 5:
  return rp;


def algorithm2(g, r_tilde, m, c = 1):

  """ @brief  Recovers r, assuming r_tilde is such that r = d * r_tilde where d
              is cm-smooth.

      This function implements Algorithm 2 from [E22p].

      [E22p] Ekerå, M.: "On the success probability of quantum order finding".
                        ArXiv 2201.07791 (2022).

      As in [E22p], d is said to be cm-smooth if d = p1^e1 * .. pk^ek, for
      q1, .., qk pairwise distinct primes, and e1, .., ek positive integer
      exponents, if it holds that qi^ei <= cm for all i in [1, k].

      @param g  The element g of order r.

      @param r_tilde  The integer r_tilde.

      @param m  A positive integer m such that r < 2^m.

      @param c  A parameter c >= 1 that specifies the maximum size of the
                missing cm-smooth component d in r = d * r_tilde.

      @return   The order r, assuming r_tilde is such that r = d * r_tilde where
                d is cm-smooth. Otherwise, None or some positive integer
                multiple of r, is returned. """

  # Step 1:
  if not is_valid_r_tilde(r_tilde, m):
    # Step 1.1:
    return None;

  # Step 2:
  x = g ** r_tilde;

  # Step 3:
  if x == 1:
    return r_tilde;

  # Step 4:
  S = [];

  # Step 5:
  for q in prime_range(floor(c * m) + 1):
    # Step 5.1:
    e = floor(log(c * m) / log(q));

    # Step 5.2:
    S.append([x, q, e]);

    # Step 5.3:
    q_pow_e = q ** e;
    x = x ** q_pow_e;

    # Step 5.4:
    if x == 1:
      # Step 5.4.1:
      break;

  # Step 6:
  if x != 1:
    # Step 6.1:
    return None;

  # Step 7:
  d = 1;

  # Step 8:
  while S != []:
    [x, q, e] = S.pop();

    # Step 8.1:
    x = x ** d;

    # Step 8.2:
    for i in range(1, e + 1):
      # Step 8.2.1:
      if x == 1:
        # Step 8.2.1.1:
        break;

      # Step 8.2.2:
      x = x ** q;
      d = d * q;

  # Step 9:
  return d * r_tilde;

# Algorithm 3, implemented as described in the paper:
def algorithm3(g, r_tilde, m, c = 1):

  """ @brief  Recovers r, assuming r_tilde is such that r = d * r_tilde where d
              is cm-smooth.

      This function implements Algorithm 3 from [E22p].

      [E22p] Ekerå, M.: "On the success probability of quantum order finding".
                        ArXiv 2201.07791 (2022).

      As in [E22p], d is said to be cm-smooth if d = p1^e1 * .. pk^ek, for
      q1, .., qk pairwise distinct primes, and e1, .., ek positive integer
      exponents, if it holds that qi^ei <= cm for all i in [1, k].

      @remark   As is explained in [E22p], this algorithm is equivalent to
                Algorithm 2 in that it computes the same result. It has a
                better worst-case time complexity than Algorithm 2, but it has
                a much higher time complexity in practice for many problem
                instances. This explains why both algorithms are included.

      @param g  The element g of order r.

      @param r_tilde  The integer r_tilde.

      @param m  A positive integer m such that r < 2^m.

      @param c  A parameter c >= 1 that specifies the maximum size of the
                missing cm-smooth component d in r = d * r_tilde.

      @return   The order r, assuming r_tilde is such that r = d * r_tilde where
                d is cm-smooth. Otherwise, None or some positive integer
                multiple of r, is returned. """

  # Step 1:
  if not is_valid_r_tilde(r_tilde, m):
    # Step 1.1:
    return None;

  # Step 2:
  def recursive(x, F):
    l = len(F);

    # Step 2.1:
    if l == 1:
      # Step 2.1.1:
      return {(F[0], x)};

    # Step 2.2:
    F_L = F[:floor(l/2)]; F_R = F[floor(l/2):];

    # Step 2.3:
    d_L = mpz(prod([q ** (floor(log(c * m) / log(q))) for q in F_R]));
    d_R = mpz(prod([q ** (floor(log(c * m) / log(q))) for q in F_L]));

    x_L = x ** d_L; x_R = x ** d_R;

    # Step 2.4:
    return recursive(x_L, F_L).union(recursive(x_R, F_R));

  # Step 3:
  x = g ** r_tilde; d = 1;

  # Step 4:
  T = recursive(x, [q for q in prime_range(floor(c * m) + 1)]);

  # Step 5:
  for (q_i, x_i) in T:
    # Step 5.1:
    e_i = 0; e_max = floor(log(c * m) / log(q_i));

    # Step 5.2:
    while x_i != 1:
      # Step 5.2.1:
      if e_i == e_max:
        return None;

      # Step 5.2.2:
      x_i = x_i ** q_i; d = d * q_i; e_i = e_i + 1;

  # Step 6:
  return d * r_tilde;

# Algorithm 4, implemented as described in the paper:
def algorithm4(g, S, m, c = 1):

  """ @brief  Returns a subset Sp of S consisting of all r_tilde in S that
              are such that d * r_tilde is a positive integer multiple of r,
              where d is cm-smooth.

      This function implements Algorithm 4 from [E22p].

      [E22p] Ekerå, M.: "On the success probability of quantum order finding".
                        ArXiv 2201.07791v2 (2022).

      As in [E22p], d is said to be cm-smooth if d = p1^e1 * .. pk^ek, for
      q1, .., qk pairwise distinct primes, and e1, .., ek positive integer
      exponents, if it holds that qi^ei <= cm for all i in [1, k].

      @param g  The element g of order r.

      @param S  A set S of candidates for the integer r_tilde.

      @param m  A positive integer m such that r < 2^m.

      @param c  A parameter c >= 1 that specifies the maximum size of the
                missing cm-smooth component d in r = d * r_tilde.

      @return   A subset Sp of S consisting of all r_tilde in S that are such
                that d * r_tilde is a positive integer multiple of r, where d
                is cm-smooth. """

  # Step 1:
  e = prod([q ** floor(log(c * m) / log(q))
    for q in prime_range(floor(c * m) + 1)]);
  x = g ** e;

  # Step 2:
  Sp = set();

  # Step 3:
  for tilde_rip in S:
    # Step 3.1:
    if is_valid_r_tilde(tilde_rip, m) and (x ** tilde_rip == 1):
      # Step 3.1.1:
      Sp.add(tilde_rip);

  # Step 4:
  return Sp;