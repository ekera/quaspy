""" @brief  A module for splitting N = pq into the large l-bit prime factors
            p and q given d = (p - 1) / 2 + (q - 1) / 2 - 2^(l - 1). """

from gmpy2 import mpz;

from gmpy2 import isqrt_rem as sqrt_rem;
from gmpy2 import ceil;

def split_N_given_d(
  d : int | mpz,
  N : int | mpz) -> set[int] | None:

  """ @brief  Splits N = pq into the large l-bit prime factors p and q given
              d = (p - 1) / 2 + (q - 1) / 2 - 2^(l - 1).

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

      This convenience function performs the last step in the above procedure.

      [EH17] Ekerå, M. and Håstad, J.: "Quantum Algorithms for Computing Short
                                        Discrete Logarithms and Factoring RSA
                                        Integers.". In: PQCrypto 2017.
                                       Springer LNCS 10346, pp. 347–363 (2017).

      [E20] Ekerå, M.: "On post-processing in the quantum algorithm for
                        computing short discrete logarithms".
                       Des. Codes Cryptogr. 88, pp. 2313–2335 (2020).

      @param d  The logarithm d.

      @param N  The integer N. It is required that N is odd, and the product
                of two large random distinct l-bit prime numbers.

      @return   A set {p, q} of two non-trivial factors of N such that N = pq,
                or None if the algorithm failed to split N. """

  # Convert inputs.
  d = mpz(d);
  N = mpz(N);

  # Sanity checks.
  if (N % 2) == 0:
    raise Exception("Error: Incorrect parameters: N is not odd.");

  if not (0 < d < N / 2):
    raise Exception("Error: Incorrect parameters: d is not on [1, N/2).");

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

  # Form d' = p + q given d.
  p_plus_q = 2 * d + 2 ** l + 2;

  # Solve d' = p + q and N = pq for p and q using the quadratic formula:
  if ((p_plus_q ** 2) - 4 * N <= 0) or (0 != (p_plus_q % 2)):
    # Failed to factor.
    return None;

  (root, remainder) = sqrt_rem((p_plus_q ** 2) - 4 * N);

  # Sanity checks.
  if (0 != remainder) or (0 != (root % 2)):
    # Failed to factor.
    return None;

  # Compute p and q.
  p = (p_plus_q - root) // 2;
  q = (p_plus_q + root) // 2;

  # Sanity checks.
  if (1 < p < N) and (1 < q < N) and (p * q == N):

    # Factored N into p and q.
    return set([int(p), int(q)]);

  # Failed to factor.
  return None;