""" @brief  A module for finding the solution to a set of congruence relations
            with coprime moduli via the Chinese remainder theorem. """

from gmpy2 import mpz;

from gmpy2 import invert;

from math import prod;

def crt(
    values : list[int | mpz],
    moduli : list[int | mpz]) -> int | mpz:

  """ @brief  Given values = [v1, ..., vn] and moduli = [N1, ..., Nn], this
              function returns an integer v in [0, N) such that v = vi (mod Ni)
              for all i in [1, n], where N = N1 * ... * Nn.

      @remark It is required that all Ni >= 2, and pairwise coprime.

      @param values   The values = [v1, ..., vn].

      @param moduli   The moduli = [N1, ..., Nn].

      @return   An integer v in [0, N) such that v = vi (mod Ni) for all i in
                [1, n], where N = N1 * ... * Nn. """

  n = len(values);

  if (n != len(moduli)) or (n < 1):
    raise Exception("Error: Incorrect parameters.");

  N = prod([Ni for Ni in moduli]);

  Tis = [N // Ni for Ni in moduli];

  Tinvis = [invert(Tis[i], moduli[i]) for i in range(n)];

  return sum([values[i] * Tinvis[i] * Tis[i] for i in range(n)]) % N;