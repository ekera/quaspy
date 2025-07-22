""" @brief  A module for computing prime ranges and prime power products, for
            sampling random primes and for testing if integers are smooth. """

from gmpy2 import mpz;

from gmpy2 import is_prime; # Imported also to be available to other modules.

from math import ceil;

from .random import sample_l_bit_integer;

def prime_range(B : int) -> list[int]:

  """ @brief  Returns an ordered list of all primes less than B.

      @param B  The upper bound B.

      @return   An ordered list of all primes less than B. """

  sieve = [True for _ in range(0, B)];
  primes = [];

  p = 2;
  while p < B:
    primes.append(p);

    for n in range(2, ceil(B / p)):
      sieve[n * p] = False;

    p += 1;
    while (p < B) and (not sieve[p]):
      p += 1;

  return primes;


def prime_power_product(B : int) -> int | mpz:

  """ @brief  Returns the product of q^e, as q runs over all primes <= B, for e
              the largest non-negative integer exponent such that q^e <= B.

      @param B  The upper bound B.

      @return   The product of q^e, as q runs over all primes <= B, for e the
                largest non-negative integer exponent such that q^e <= B."""

  factor = 1;

  for q in prime_range(B + 1):
    q_pow_e = q;

    while True:
      tmp = q_pow_e * q;
      if tmp > B:
        break;
      q_pow_e = tmp;

    factor *= q_pow_e;

  return factor;


def is_B_smooth(d : int | mpz, B : int) -> bool:

  """ @brief  Tests if the integer d is B-smooth.

      As in [E24], d is said to be B-smooth if d = p1^e1 * .. pk^ek, for
      q1, ..., qk pairwise distinct primes, and e1, ..., ek positive integer
      exponents, if it holds that qi^ei <= B for all i in [1, k].

      [E24] Ekerå, M.: "On the success probability of quantum order finding".
                       ACM Trans. Quantum Comput. 5(2):11 (2024).

      @param d  The integer d.

      @param B  The upper bound B.

      @return   True if d is B-smooth, False otherwise. """

  if d <= 0:
    raise Exception("Error: Incorrect parameters.");

  for q in prime_range(B + 1):
    q_pow_e = 1;

    while (d % q) == 0:
      q_pow_e *= q;
      if q_pow_e > B:
        return False;
      d //= q;

  return d == 1;


def sample_l_bit_prime(l : int) -> int | mpz:

  """ @brief  Returns an l-bit prime selected uniformly at random from the set
              of all such primes.

      @remark   This function calls sample_l_bit_integer() to select l-bit
                integers uniformly at random from the set of all l-bit integers.
                In practice, sample_l_bit_integer() returns an l-bit integer
                that may be conjectured to be indistinguishable from an integer
                that is selected uniformly at random from said set.

      @remark   This function assumes that the is_prime() function performs a
                deterministic primality test. In practice, this test is likely
                probabilistic, but the probability of incorrectly classifying
                an composite as a prime is so small that the test may be
                conjectured to be indistinguishable from a determinstic test.

      @param l  The bit length l of the prime to sample.

      @return   An l-bit prime selected uniformly at random from the set of all
                such primes. """

  while True:
    p = mpz(sample_l_bit_integer(l));
    if p.is_prime():
      return p;