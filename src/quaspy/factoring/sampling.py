""" @brief  A module for exactly or heuristically sampling an element g
            uniformly at random from the multiplicative group of the ring of
            integers modulo N and returning the order r of g. This when the
            complete factorization of N is known.

    The procedures in this module are described in [E21b], in [E24t] (see in
    particular Sect. 5.2.3) and in the Factoritall repository (available at
    https://github.com/ekera/factoritall).

    [E21b] Ekerå, M.: "On completely factoring any integer efficiently in a
                       single run of an order-finding algorithm".
                      Quantum Inf. Process. 20(6):205 (2021).

    [E24t] Ekerå, M.: "On factoring integers, and computing discrete logarithms
                       and orders, quantumly".
                      PhD thesis, KTH Royal Institute of Technology (2024). """

from gmpy2 import mpz;
from gmpy2 import powmod;

from gmpy2 import gcd;
from gmpy2 import lcm;

from math import prod;
from math import floor;

from ..math.crt import crt;
from ..math.primes import prime_range;
from ..math.random import sample_integer;

B_DEFAULT_SAMPLE = 10 ** 6;

def sample_r_given_N(
  N : int | mpz,
  factors : list[list[int | mpz]]) -> int:

  """ @brief  Returns the order r of an element g selected uniformly at random
              from the multiplicative group of the ring of integers modulo N,
              without explicitly computing and returning the element g.

      Suppose that N = p1^e1 * ... * pn^en, for p1, ..., pn pairwise distinct
      odd prime factors, and e1, ..., en positive integer exponents.

      For i in [1, n], suppose that gi is selected uniformly at random from the
      multiplicative group of the ring of integers modulo pi^ei, and that the
      order ri of gi is computed. Then r = lcm(r1, ..., rn) is the order of g,
      where g may be computed via the Chinese remainder theorem, by using that
      it must hold that gi = g mod pi^ei for i in [1, n].

      (Note that the above map is an isomorphism: Hence g will be selected
      uniformly at random from the multiplicative group of the ring of integers
      modulo N, as the gi are selected uniformly at random from the
      multiplicative group of the ring of integers modulo pi^ei.)

      A problem with the above approach is that it is hard to compute the order
      ri, and hence to compute r, unless the factorization of pi - 1 for i in
      [1, n] is known. To circumvent this problem, instead of directly selecting
      the gi uniformly at random from the multiplicative group of the ring of
      integers modulo pi^ei, suppose exponents di are instead selected uniformly
      at random on the interval [0, lambda(pi^ei)), where

          lambda(pi^ei) = (pi - 1) pi^(ei - 1)

      as all pi are odd, and gi = Gi^di computed, for Gi some fixed generator of
      the multiplicative group the ring of integers modulo pi^ei. Then, gi is of
      order ri = lambda(pi^ei) / gcd(lambda(pi^ei), di) for i in [1, n], so the
      ri are easy to compute, as is r = lcm(r1, ..., rn).

      Again, unless the factorization of pi - 1 for i in [1, n] is known, it is
      hard to prove that an element Gi is a generator, and hence to explicitly
      compute g. This explains why this function only returns r.

      The above procedure is not described in [E21b], but in [E24t] (see in
      particular Sect. 5.2.3) and in the Factoritall repository (available at
      https://github.com/ekera/factoritall).

      [E21b] Ekerå, M.: "On completely factoring any integer efficiently in a
                         single run of an order-finding algorithm".
                        Quantum Inf. Process. 20(6):205 (2021).

      [E24t] Ekerå, M.: "On factoring integers, and computing discrete
                         logarithms and orders, quantumly". PhD thesis, KTH
                        Royal Institute of Technology (2024).

      @param N  The integer N.

      @param factors  The factors of N = p1^e1 * ... * pn^en, represented on the
                      form [[p1, e1], ..., [pn, en]], for p1, ..., pn pairwise
                      distinct prime factors, and for e1, ..., en positive
                      integer exponents.

      @return   The order r of an element g selected uniformly at random from
                the multiplicative group of the ring of integers modulo N. """

  # Make sure the factors are GMP integers.
  factors = [[mpz(pi), int(ei)] for [pi, ei] in factors];
  n = len(factors);

  # Sanity checks.
  if N != prod([pi ** ei for [pi, ei] in factors]):
    raise Exception("Error: Incorrect factorization: N != sum(pi ** ei)");

  pis = [pi for [pi, _] in factors];
  if len(set(pis)) != n:
    raise Exception("Error: Incorrect factorization: " +
      "The pis are not all pairwise distinct.");

  if min(pis) <= 2:
    raise Exception("Error: Incorrect factorization: " +
      "It is required that all pi > 2.");

  for pi in pis:
    if not pi.is_prime():
      raise Exception("Error: Incorrect factorization: " +
        "It is required that all pis be prime.");

  eis = [ei for [_, ei] in factors];
  if min(eis) < 1:
    raise Exception("Error: Incorrect factorization: " +
      "It is required that all ei >= 1.");

  # Compute lambdai for i in [1, n].
  lambdais = [(pi - 1) * (pi ** (ei - 1)) for [pi, ei] in factors];

  # Sample di for i in [1, n].
  dis = [sample_integer(lambdai) for lambdai in lambdais];

  # Compute ri for i in [1, n].
  ris = [lambdais[i] // gcd(lambdais[i], dis[i]) for i in range(n)];

  # Compute r.
  r = ris[0];
  for i in range(1, n):
    r = lcm(r, ris[i]);

  # Return r.
  return int(r);


def sample_g_r_given_N(
  N : int | mpz,
  N_factors : list[list[int | mpz]],
  pi_minus_one_factors : list[list[list[int | mpz]]] | None = None,
  B : int = B_DEFAULT_SAMPLE) -> list[int]:

  """ @brief  Returns [g, r], for g an element selected uniformly at random
              from the multiplicative group of the ring of integers modulo N,
              and r either a heuristic estimate of the order of g, or the exact
              order of g, depending on if optional parameters are specified.

      Suppose that N = p1^e1 * ... * pn^en, for p1, ..., pn pairwise distinct
      odd prime factors, and e1, ..., en positive integer exponents.

      For i in [1, n], this function then selects gi from the multiplicative
      group of the ring of integers modulo pi^ei.

      1. If the factorization of pi - 1 for i in [1, n] is *not* specified, this
      function then heuristically estimates the order ri of gi by using the
      method in App. A of [E21b]: Specificially, by using that

         lambda(pi^ei) = (pi - 1) pi^(ei - 1),

      as pi is odd, and by using a factor base of primes <= B to find all small
      factors of pi - 1 via trial division. It then computes g via the Chinese
      remainder theorem, by requiring that gi = g mod pi^ei, along with a
      heuristic estimate r = lcm(r1, ..., rn) of the order of g, that is
      correct with high probability, as is explained in [E21b].

      2. If the factorization of pi - 1 for i in [1, n] is specified, this
      function then exactly computes the order ri of gi. Specifically, by using

          lambda(pi^ei) = (pi - 1) pi^(ei - 1)

      as an initial guess ri' for the order ri of gi. Then, for each prime
      factor f that divide ri', for as long as f divides ri' and
      gi^(ri' / f) = 1 (mod N), let ri' <- ri' / f. It follows that ri = ri' at
      the end of the procedure. The order of g is then r = lcm(r1, ..., rn).

      The above procedure is described in [E21b], and in the factoritall
      repository (available at https://github.com/ekera/factoritall).

      [E21b] Ekerå, M.: "On completely factoring any integer efficiently in a
                         single run of an order-finding algorithm".
                         Quantum Inf. Process. 20(6):205 (2021).

      @param N  The integer N.

      @param N_factors  The factors of N = p1^e1 * ... * pn^en, represented on
                        the form [[p1, e1], ..., [pn, en]], for p1, ..., pn
                        pairwise distinct prime factors, and for e1, ..., en
                        positive integer exponents.

      @param pi_minus_one_factors  The factors of pi-1 = qi1^di1 * ... * qim^dim,
                                   for i in [1, n], represented on the form
                                   [F1, ..., Fn], where each Fi is on the form
                                   [[qi1, qi1], ..., [qim, qim]], for
                                   qi1, ..., qim pairwise distinct prime factors,
                                   and for di1, ..., dim positive integer
                                   exponents. May be set to None, in which case
                                   r will be computed deterministically as
                                   described above. If explicitly specified, the
                                   order r will be computed exactly.

      @param B  The upper bound on the prime factors to consider when performing
                trial division. Has no effect if pi_minus_one_factors is
                explicitly specified, as trial division is then not performed.

      @return   The pair [g, r], for g an element selected uniformly at random
                from the multiplicative group of the ring of integers modulo N,
                and r a heuristic estimate of the order of g, or the exact order
                of g, depending on if optional parameters are specified. """

  # Process the factors [pi, ei] of N:

  # Make sure the factors are GMP integers.
  factors = [[mpz(pi), int(ei)] for [pi, ei] in N_factors];

  n = len(factors);

  # Sanity checks.
  if N != prod([pi ** ei for [pi, ei] in factors]):
    raise Exception("Error: Incorrect factorization of N: N != sum(pi ** ei)");

  pis = [pi for [pi, _] in factors];
  if len(set(pis)) != n:
    raise Exception("Error: Incorrect factorization of N: " +
      "The pis are not all pairwise distinct.");

  if min(pis) <= 2:
    raise Exception("Error: Incorrect factorization of N: " +
      "It is required that all pi > 2.");

  for pi in pis:
    if not pi.is_prime():
      raise Exception("Error: Incorrect factorization of N: " +
        "It is required that all pis be prime.");

  eis = [ei for [_, ei] in factors];
  if min(eis) < 1:
    raise Exception("Error: Incorrect factorization of N: " +
      "It is required that all ei >= 1.");

  if pi_minus_one_factors != None:
    # Also process the factors [qj, dj] of pi - 1:

    if len(pi_minus_one_factors) != n:
      raise Exception("Error: Incorrect factorization of pi - 1: "
        "Expected one entry for each of the n prime factors pi of N.");

    for i in range(n):
      # Make sure the factors are GMP integers.
      pi_minus_one_factors[i] = \
        [[mpz(qj), dj] for [qj, dj] in pi_minus_one_factors[i]];

      mi = len(pi_minus_one_factors[i]);

      # Sanity checks.
      if pis[i] - 1 != prod([qj ** dj for [qj, dj] in pi_minus_one_factors[i]]):
        raise Exception("Error: Incorrect factorization of pi - 1: "\
          "p" + str(i) + " - 1 != sum(qj ** dj)");

      qjs = [qj for [qj, _] in pi_minus_one_factors[i]];
      if len(set(qjs)) != mi:
        raise Exception("Error: Incorrect factorization of pi - 1: " +
          "The qjs are not all pairwise distinct.");

      if min(qjs) < 2:
        raise Exception("Error: Incorrect factorization of pi - 1: " +
          "It is required that all qj >= 2.");

      for qj in qjs:
        if not qj.is_prime():
          raise Exception("Error: Incorrect factorization of pi - 1: " +
            "It is required that all qjs be prime.");

      djs = [dj for [_, dj] in pi_minus_one_factors[i]];
      if min(djs) < 1:
        raise Exception("Error: Incorrect factorization of pi - 1: " +
          "It is required that all dj >= 1.");

  # Pre-compute the set of all primes on [2, B).
  primes = prime_range(B);

  # Heuristically compute gi and ri for i in [1, n].
  gis = [];
  ris = [];

  # Compute modulii = pi^ei for i in [1, n].
  moduliis = [pi ** ei for [pi, ei] in factors];

  for i in range(n):
    pi = pis[i];
    ei = eis[i];
    modulii = moduliis[i];

    # Sample gi in Z_{pi^ei}^*.
    while True:
      gi = sample_integer(modulii);
      if gcd(gi, pi) == 1:
        break;

    # Factor pi - 1.
    if pi_minus_one_factors != None:
      F = pi_minus_one_factors[i];
      if ei > 1:
        F.append([pi, ei - 1]);

      ri = 1;
      gip = gi;
    else:
      F = [];

      if ei > 1:
        F.append([pi, ei - 1]);

      tmp = pi - 1;

      for q in primes:
        d = 0;

        while (tmp % q) == 0:
          d += 1;
          tmp //= q;

        if d > 0:
          F.append([q, d]);

      ri = tmp;
      gip = powmod(gi, tmp, N);

    def recursive(x, F):
      l = len(F);

      if l == 0:
        return set();

      if l == 1:
        [q, d] = F[0];
        return {(x, q, d)};

      F_L = F[:floor(l/2)]; F_R = F[floor(l/2):];

      d_L = mpz(prod([q ** d for [q, d] in F_R]));
      d_R = mpz(prod([q ** d for [q, d] in F_L]));

      x_L = powmod(x, d_L, modulii); x_R = powmod(x, d_R, modulii);

      return recursive(x_L, F_L).union(recursive(x_R, F_R));

    for (x, q, d) in recursive(gip, F):
      for i in range(d):
        if x == 1:
          break;

        x = powmod(x, q, modulii);
        ri *= q;

      if x != 1:
        raise Exception("Error: Sanity check failed: Internal error.")

    gis.append(gi);
    ris.append(ri);

  # Compute r.
  r = ris[0];
  for i in range(1, n):
    r = lcm(r, ris[i]);

  # Compute g.
  g = crt(gis, moduliis);

  # Sanity checks.
  for i in range(n):
    if gis[i] != (g % moduliis[i]):
      raise Exception("Error: Sanity check failed: Internal error.");

  if powmod(g, r, N) != 1:
    raise Exception("Error: Sanity check failed: Internal error.");

  for q in set(primes + pis):
    if (r % q) == 0:
      if powmod(g, r // q, N) == 1:
        raise Exception("Error: Sanity check failed: Internal error.");

  # Return [g, r].
  return [int(g), int(r)];