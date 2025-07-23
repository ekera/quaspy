""" @brief  A module for splitting N, given the order r of an element g selected
            uniformly at random from the multiplicative group of the ring of
            integers modulo N, where g must be explicitly specified. This by
            using the original algorithm from [Shor94].

    [Shor94] Shor, P.W.: "Algorithms for Quantum Computation: Discrete
                          Logarithms and Factoring".
                         In: Proceedings from FOCS '94, pp. 124–134 (1994). """

from gmpy2 import mpz;

from gmpy2 import gcd;
from gmpy2 import powmod;

def split_N_given_g_r(
  g : int | mpz,
  r : int | mpz,
  N : int | mpz) -> set[int] | None:

  """ @brief  Attempts to split N, given the order r of an element g selected
              uniformly at random from the multiplicative group of the ring of
              integers modulo N, where g must be explicitly specified. This by
              using the original algorithm from [Shor94].

      [Shor94] Shor, P.W.: "Algorithms for Quantum Computation: Discrete
                            Logarithms and Factoring".
                           In: Proceedings from FOCS '94, pp. 124–134 (1994).

      The algorithm succeeds iff r is even and g^(r/2) != -1 (mod N).

      @param g  The element g.

      @param r  The order r of g mod N.

      @param N  The modulus N. As in [Shor94], it is required that N is odd,
                and not a perfect prime power.

      @return   A set {p, q} of two non-trivial factors of N such that N = pq,
                or None if the algorithm failed to split N. """

  # Convert inputs.
  g = mpz(g);
  r = mpz(r);
  N = mpz(N);

  # Sanity checks.
  if not (0 < g < N):
    raise Exception("Error: Incorrect parameters: g is not on [1, N).");

  if (N % 2) == 0:
    raise Exception("Error: Incorrect parameters: N is not odd.");

  if not (0 < r < N / 2):
    raise Exception("Error: Incorrect parameters: N is not on [1, N/2).");

  if not gcd(g, N) == 1:
    raise Exception("Error: Incorrect parameters: g is not coprime to N.");

  # Attempt to split N.
  if (r % 2) != 0:
    # Additional sanity check.
    if powmod(g, r, N) != 1:
      raise Exception("Error: Incorrect parameters: r is not the order of g.");

    return None;

  x = powmod(g, r // 2, N);

  # Additional sanity check.
  if (x == 1) or (powmod(x, 2, N) != 1):
    raise Exception("Error: Incorrect parameters: r is not the order of g.");

  if x == N - 1:
    return None;

  [p, q] = [gcd(x - 1, N), gcd(x + 1, N)];

  # Additional sanity check.
  if N != p * q:
    raise Exception("Error: Internal error: Is N a perfect prime power?");

  return set([int(p), int(q)]);