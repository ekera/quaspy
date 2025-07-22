""" @brief  A module for factoring large random RSA integers.

    This module uses Ekerå–Håstad's algorithm to factor RSA integers, as
    described in [EH17], with improvements from [E20] and [E23p].

      [EH17] Ekerå, M. and Håstad, J.: "Quantum Algorithms for Computing Short
                                        Discrete Logarithms and Factoring RSA
                                        Integers.". In: PQCrypto 2017.
                                       Springer LNCS 10346, pp. 347–363 (2017).

      [E20] Ekerå, M.: "On post-processing in the quantum algorithm for
                        computing short discrete logarithms".
                       Des. Codes Cryptogr. 88, pp. 2313–2335 (2020).

      [E23p] Ekerå, M.: "On the success probability of the quantum algorithm for
                         the short DLP". ArXiv 2309.01754v2 (2025). """

from gmpy2 import mpz;

from gmpy2 import ceil;

from ...math.groups import CyclicGroupElement;

def setup_x_given_g_N(
    g : CyclicGroupElement,
    N : int | mpz) -> CyclicGroupElement:

  """ @brief  Sets up x = g^d' for d' = (N - 1) / 2 - 2^(l - 1) given g and N,
              for N the product of two large random distinct l-bit primes.

      This is a convenience function for Ekerå–Håstad's algorithm [EH17] that
      factors a large random RSA integer N = pq into p and q.

      More specifically, as is explained in App. A.2 of [E20], it holds that

        x = g^d' = g^((N - 1) / 2 - 2^(l - 1))
                 = g^((p - 1) / 2 + (q - 1) / 2 - 2^(l - 1)) = g^d.

      Provided that the order r of g is sufficiently large, it holds that

        d = d' mod r = (p - 1) / 2 + (q - 1) / 2 - 2^(l - 1),

      allowing d to be computed as the discrete logarithm of x to the base g.

      This may be done efficiently by using Ekerå–Håstad's quantum algorithm
      that computes short discrete logarithms in groups of unknown order.

      Given d = (p - 1) / 2 + (q - 1) / 2 - 2^(l - 1) and N = pq, it is then
      trivial to compute p and q by solving a quadratic equation.

      This convenience function computes x in the above procedure given g and N.

      [EH17] Ekerå, M. and Håstad, J.: "Quantum Algorithms for Computing Short
                                        Discrete Logarithms and Factoring RSA
                                        Integers.". In: PQCrypto 2017.
                                       Springer LNCS 10346, pp. 347–363 (2017).

      [E20] Ekerå, M.: "On post-processing in the quantum algorithm for
                        computing short discrete logarithms".
                       Des. Codes Cryptogr. 88, pp. 2313–2335 (2020).

      @param g  The generator g selected uniformly at random from the
                multiplicative group of the ring of integers modulo N.

      @param N  The integer N = pq, for p and q two large distinct random
                l-bit prime numbers.

      @return   The element x = g^((N - 1) / 2 - 2^(l - 1)). """

  # Convert inputs.
  N = mpz(N);

  # Sanity checks.
  if 1 != (N % 2):
    raise Exception("Error: Expected N to be odd.");

  # Compute the bit length l of p and q where N = pq: We have that
  #
  #   p, q >= 2^(l-1) => N = pq = 2^(2(l-1))
  #                        = 2^(2l-2),
  #
  # i.e. N is an n >= 2l - 1 bit integer. At the same time, we have that
  #
  #   p, q < 2^l => N = pq < 2^2l,
  #
  # i.e. N is an n <= 2l bit integer. Hence, for n the bit length of N, it must
  # be that l = ceil(n / 2), so we may compute l from n in this manner:

  l = int(ceil(N.bit_length() / 2));

  # Select the exponent d' as explained in [E20] and compute x.
  x = g ** ((N - 1) // 2 - 2 ** (l - 1));

  # Return x.
  return x;


def setup_d_given_p_q(
  p : int | mpz,
  q : int | mpz) -> int:

  """ @brief  Sets up d = (p - 1) / 2 + (q - 1) / 2 - 2^(l - 1) given p and q.

      This is a convenience function for Ekerå–Håstad's algorithm [EH17] that
      factors a large random RSA integer N = pq into p and q.

      More specifically, as is explained in App. A.2 of [E20], it holds that

        x = g^d' = g^((N - 1) / 2 - 2^(l - 1))
                 = g^((p - 1) / 2 + (q - 1) / 2 - 2^(l - 1)) = g^d

      Provided that the order r of g is sufficiently large, it holds that

        d = d' mod r = (p - 1) / 2 + (q - 1) / 2 - 2^(l - 1),

      allowing d to be computed as the discrete logarithm of x to the base g.

      This may be done efficiently by using Ekerå–Håstad's quantum algorithm
      that computes short discrete logarithms in groups of unknown order.

      Given d = (p - 1) / 2 + (q - 1) / 2 - 2^(l - 1) and N = pq, it is then
      trivial to compute p and q by solving a quadratic equation.

      This convenience function computes d in the above procedure given p and q.

      [EH17] Ekerå, M. and Håstad, J.: "Quantum Algorithms for Computing Short
                                        Discrete Logarithms and Factoring RSA
                                        Integers.". In: PQCrypto 2017.
                                       Springer LNCS 10346, pp. 347–363 (2017).

      [E20] Ekerå, M.: "On post-processing in the quantum algorithm for
                        computing short discrete logarithms".
                       Des. Codes Cryptogr. 88, pp. 2313–2335 (2020).

      @param p  A large random l-bit prime.
      @param q  A large random l-bit prime not equal to p.

      @return   The logarithm d = ((p - 1) / 2) + ((q - 1) / 2) - 2^(l - 1). """

  # Sanity checks.
  if p == q:
    raise Exception("Error: The primes p and q must be distinct.");

  if ((p % 2) != 1) or ((q % 2) != 1):
    raise Exception("Error: The primes p and q must be odd.");

  # Compute the bit length l of p and q.
  l = p.bit_length();
  if q.bit_length() != l:
    raise Exception("Error: The primes p and q must be of equal bit length.");

  # Compute d as explained in [E20].
  d = ((p - 1) // 2) + ((q - 1) // 2) - 2 ** (l - 1);

  # Return d;
  return int(d);