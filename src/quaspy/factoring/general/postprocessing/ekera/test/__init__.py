""" @brief  A module for testing the functions in the parent module for
            factoring N completely.

    This module implements the test procedure described in App. A to [E21b].

    [E21b] Ekerå, M.: "On completely factoring any integer efficiently in a
                       single run of an order-finding algorithm".
                      Quantum Inf. Process. 20(6):205 (2021).

    @remark This module is modelled upon the corresponding Sage script in the
            Factoritall repository, see https://github.com/ekera/factoritall.

    To run all unit tests in this module by calling test_all(), execute:

    ```sh

    $ python3 -m quaspy.factoring.general.postprocessing.ekera.test

    ``` """

from gmpy2 import mpz;

from math import gcd;
from math import prod;
from math import floor;
from math import prod;

from ......utils.timer import Timer;

from ......math.groups import SimulatedCyclicGroupElement;
from ......math.groups import IntegerModRingMulSubgroupElement;

from ......math.random import sample_integer;

from .. import solve_j_for_factors;
from .. import solve_r_for_factors;

from ......orderfinding.general.sampling import sample_j_given_r;

from .....sampling import sample_r_given_N;
from .....sampling import sample_g_r_given_N;

from .....sampling import B_DEFAULT_SAMPLE;

from ......orderfinding.general.postprocessing.ekera import B_DEFAULT_SOLVE;

from ......math.primes import is_B_smooth;

def is_admissible_r(
  r : int | mpz,
  B : int,
  N_factors : list[list[int | mpz]]) -> bool:

  """ @brief  Tests if (pi - 1) / gcd(pi - 1, r) is B-smooth for all but at
              most one of the prime factors pi of N.

      @param r  The order r of an

      @param B  The upper bound B.

      @param N_factors  The factors of N = p1^e1 * ... * pn^en, represented on
                        the form [[p1, e1], ..., [pn, en]], for p1, ..., pn
                        pairwise distinct prime factors, and for e1, ..., en
                        positive integer exponents.

      @return   True if (pi - 1) / gcd(pi - 1, r) is B-smooth for all but at
                most one of the prime factors pi of N, False otherwise. """

  B = floor(B);

  count = 0;

  for pi in [pi for [pi, _] in N_factors]:
    # Compute the missing factor for each prime.
    di = (pi - 1) // gcd(pi - 1, r);

    # Check if the missing factor is B-smooth. If not, add to count.
    if not is_B_smooth(di, B):
      count += 1;

  return count < 2;


def test_solve_r_for_factors_of_N(
  N : int | mpz,
  N_factors : list[list[int | mpz]],
  c : int = 1,
  verbose : bool = True) -> None:

  """ @brief  Samples the order r of an element g selected uniformly at random
              from the multiplicative group of the ring of integers modulo N,
              and attempts to recover the complete factorization of N given r.

      More specifically, this function accepts as input the complete
      factorization of an integer N = p1^e1 * ... * pn^en with n distinct prime
      factors pi, and optionally of pi - 1, for i in [1, n].

      It then samples the order r of an element g selected uniformly at random
      from the multiplicative group of the ring of integers modulo N, and
      attempts to solve r for the complete factorization of N.

      @param N  The integer N.

      @param N_factors  The factors of N = p1^e1 * ... * pn^en, represented on
                        the form [[p1, e1], ..., [pn, en]], for p1, ..., pn
                        pairwise distinct prime factors, and for e1, ..., en
                        positive integer exponents.

      @param c  A parameter c >= 1 that specifies the maximum size of the
                missing smooth component in lambda'(N) when solving r for the
                complete factorization of N.

      @param verbose  A flag that may be set to True to print intermediary
                      results. Defaults to True.

      @return   This function has no return value. If the test fails, or if some
                other error occurs, an exception is instead raised. """

  # Sample r.
  while True:
    r = sample_r_given_N(N, N_factors);
    if r != 1:
      break;

  if verbose:
    print("Sampled r = " + str(r) + "\n");

  # Test r.
  if not is_admissible_r(r, c * N.bit_length(), N_factors):
    if verbose:
      print("\n*** Sampled r which does not admit factorization:\n", d);

    return;

  if verbose:
    print("Done setting up the problem instance. Solving commences...\n");

  # Setup and start a timer.
  timer = Timer().start();

  # Attempt to recover the factors.
  candidate_factors = solve_r_for_factors(
                        r = r,
                        N = N,
                        c = c,
                        verbose = verbose);

  # Stop the timer.
  timer.stop();

  # Compare.
  if candidate_factors != set([pi for [pi, _] in N_factors]):
    raise Exception("Error: An incorrect complete factorization was returned.");

  if verbose:
    print("\n*** Successfully recovered the complete factorization.");


def test_solve_r_for_factors_of_random_N(
  l : int = 1024,
  n : int = 2,
  e_max : int = 1,
  c : int = 1,
  verbose : bool = True) -> None:

  """ @brief  Samples the order r of an element g selected uniformly at random
              from the multiplicative group of the ring of integers modulo N,
              and attempts to recover the complete factorization of N given r.

      More specifically, this function first selects n pairwise distinct l-bit
      prime factors pi from the set of all l-bit prime factors, and n integer
      exponents ei uniformly at random from the interval [1, e_max), for i on
      [1, n], after which it computes the integer N = p1^e1 * ... * pn^en.

      It then samples the order r of an element g selected uniformly at random
      from the multiplicative group of the ring of integers modulo N, and
      attempts to solve r for the complete factorization of N.

      @remark   This convenience function simply sets up the factors of N, and
                then calls the test_solve_r_for_factors_of_N() function.

      @param l  The length in bits of each distinct prime factors of N.

      @param n  The number of distinct prime factors of N.

      @param e_max  An upper bound on the exponent of each distinct prime
                    factor.

      @param c  A parameter c >= 1 that specifies the maximum size of the
                missing smooth component in lambda'(N) when solving r for the
                complete factorization of N.

      @param verbose  A flag that may be set to True to print intermediary
                      results. Defaults to True.

      @return   This function has no return value. If the test fails, or if some
                other error occurs, an exception is instead raised. """

  # Select factors.
  if verbose:
    print("Sampling N, please wait, this may take a moment...\n");

  factors = [];
  pis = set();

  for i in range(n):
    if e_max == 1:
      ei = 1;
    else:
      ei = int(1 + sample_integer(e_max - 1));

    while True:
      pi = mpz(2 ** (l - 1) + sample_integer(2 ** (l - 1)));

      if not pi.is_prime():
        continue;

      if pi in pis:
        continue;

      break;

    if verbose:
      print("Selected factor " + str(i) + ":", str(pi) + "^" + str(ei));

    factors.append([pi, ei]);

  # Compute N.
  N = prod([pi ** ei for [pi, ei] in factors]);

  if verbose:
    print("\nSelected N = " + str(N) + "\n");

  # Attempt to solve.
  test_solve_r_for_factors_of_N(
    N = N,
    N_factors = factors,
    c = c,
    verbose = verbose);


def test_solve_j_for_factors_of_N(
  N : int | mpz,
  N_factors : list[list[int | mpz]],
  pi_minus_one_factors : list[list[list[int | mpz]]] | None = None,
  c_factor : int = 1,
  c_solve : int = 1,
  B_solve : int = B_DEFAULT_SOLVE,
  B_sample : int = B_DEFAULT_SAMPLE,
  sample_g : bool = False,
  verbose : bool = True) -> None:

  """ @brief  Samples a frequency j from the probability distribution induced
              by Shor's order-finding algorithm when factoring N, and then
              attempts to solve j for the complete factorization of N.

      More specifically, this function accepts as input the complete
      factorization of an integer N = p1^e1 * ... * pn^en with n distinct prime
      factors pi, and optionally of pi - 1, for i in [1, n].

      It samples the order r of an element g selected uniformly at random from
      the multiplicative group of the ring of integers modulo N.

      Next, it samples a frequency j from the probability distribution induced
      by the quantum order-finding algorithm when performing order-finding with
      respect to an element of order r.

      Finally, it attempts to solve j for the complete factorization of N, by
      first solving j for r, and by then solving r for the complete
      factorization of N.

      @param N  The integer N.

      @param N_factors  The factors of N = p1^e1 * ... * pn^en, represented on
                        the form [[p1, e1], ..., [pn, en]], for p1, ..., pn
                        pairwise distinct prime factors, and for e1, ..., en
                        positive integer exponents.

      @param pi_minus_one_factors  The factors of

                                      pi - 1 = qi1^di1 * ... * qim^dim,

                                   for i in [1, n], represented on the form
                                   [F1, ..., Fn], where each Fi is on the form
                                   [[qi1, qi1], ..., [qim, qim]], for
                                   qi1, ..., qim pairwise distinct prime
                                   factors, and for di1, ..., dim positive
                                   integer exponents.

                                   May be set to None, in which case r will be
                                   computed deterministically as described
                                   above. If explicitly specified, the order r
                                   will be computed exactly.

      @param c_factor   A parameter c_factor >= 1 that specifies the maximum
                        size of the missing smooth component in lambda'(N) when
                        solving r for the complete factorization of N.

      @param c_solve  A parameter c_solve >= 1 that specifies the maximum size
                      of the missing smooth component d in r = d * r_tilde when
                      solving j for r, or a multiple of r.

      @param B_solve  The bound B on the offset in j when solving j for r.

      @param B_sample   The bound B on the offset in j when sampling j from the
                        frequency distribution induced by Shor's algorithm.

      @param sample_g   A flag that may be set to True to not exactly sample the
                        order r of an element g selected uniformly at random
                        from the multiplicative group of the ring of integers
                        modulo N, but rather to sample g and to then exactly or
                        heuristically compute r (depending on whether the
                        complete factorizations of pi - 1 for i in [1, n] are
                        specified or not).

      @param verbose  A flag that may be set to True to print intermediary
                      results. Defaults to True.

      @return   This function has no return value. If the test fails, or if some
                other error occurs, an exception is instead raised. """

  m = l = N.bit_length() - 1;
  while (l > 1) and (((N / 2) ** 2) < (2 ** (m + l - 1))):
    l -= 1;

  if sample_g:
    # Sample r and g.
    if verbose:
      print("Sampling g, please wait, this may take a moment...\n");

    while True:
      [g, r] = sample_g_r_given_N(N, N_factors, pi_minus_one_factors);
      if r != 1:
        break;

    if verbose:
      print("Sampled g = " + str(g) + "\n");
      print("Heuristically computed r = " + str(r) + "\n");

    # Setup g.
    g = IntegerModRingMulSubgroupElement(g, N);
  else:
    if pi_minus_one_factors != None:
      raise Exception(
        "Error: Specify pi_minus_one_factors only when sample_g = True.");

    # Sample r.
    while True:
      r = sample_r_given_N(N, N_factors);
      if r != 1:
        break;

    if verbose:
      print("Sampled r = " + str(r) + "\n");

    # Setup g.
    g = SimulatedCyclicGroupElement(r);

  # Test r.
  if not is_admissible_r(r, c_factor * N.bit_length(), N_factors):
    if verbose:
      print("\n*** Sampled r which does not admit factorization:\n", d);

    return;

  # Sample j.
  j = sample_j_given_r(r = r, m = m, l = l, B = B_sample);
  if None == j:
    if verbose:
      print("\n*** Failed to sample j.\n");

    return;

  if verbose:
    print("Sampled j = " + str(j) + "\n");

    print("Done setting up the problem instance. Solving commences...\n");

  # Setup and start a timer.
  timer = Timer().start();

  # Attempt to recover the factors.
  candidate_factors = solve_j_for_factors(
                        j = j,
                        m = m,
                        l = l,
                        g = g,
                        N = N,
                        c_factor = c_factor,
                        c_solve = c_solve,
                        B = B_solve,
                        accept_multiple = True,
                        verbose = verbose,
                        opt_speculative = True);

  # Stop the timer.
  timer.stop();

  # Compare.
  if candidate_factors != set([pi for [pi, _] in N_factors]):
    raise Exception("Error: An incorrect complete factorization was returned.");

  if verbose:
    print("\n*** Successfully recovered the complete factorization.");


def test_solve_j_for_factors_of_random_N(
  l : int = 1024,
  n : int = 2,
  e_max : int = 1,
  c_factor : int = 1,
  c_solve : int = 1,
  B_solve : int = B_DEFAULT_SOLVE,
  B_sample : int = B_DEFAULT_SAMPLE,
  sample_g : bool = False,
  verbose : bool = True) -> None:

  """ @brief  Samples a frequency j from the probability distribution induced
              by Shor's order-finding algorithm when factoring N, and then
              attempts to solve j for the complete factorization of N.

      More specifically, this function first selects n pairwise distinct l-bit
      prime factors pi from the set of all l-bit prime factors, and n integer
      exponents ei uniformly at random from the interval [1, e_max), for i on
      [1, n], after which it computes the integer N = p1^e1 * ... * pn^en.

      It then samples the order r of an element g selected uniformly at random
      from the multiplicative group of the ring of integers modulo N. Next, it
      samples a frequency j from the probability distribution induced by the
      quantum order-finding algorithm when performing order-finding with respect
      to an element of order r.

      Finally, it attempts to solve j for the complete factorization of N, by
      first solving j for r, and by then solving r for the complete
      factorization of N.

      @remark   This convenience function simply sets up the factors of N, and
                then calls the test_solve_j_for_factors_of_N() function.

      @param l  The length in bits of each distinct prime factors of N.

      @param n  The number of distinct prime factors of N.

      @param e_max  An upper bound on the exponent of each distinct prime
                    factor.

      @param c_factor   A parameter c_factor >= 1 that specifies the maximum
                        size of the missing smooth component in lambda'(N) when
                        solving r for the complete factorization of N.

      @param c_solve  A parameter c_solve >= 1 that specifies the maximum size
                      of the missing smooth component d in r = d * r_tilde when
                      solving j for r, or a multiple of r.

      @param B_solve  The bound B on the offset in j when solving j for r.

      @param B_sample   The bound B on the offset in j when sampling j from the
                        frequency distribution induced by the quantum
                        order-finding algorithm.

      @param sample_g   A flag that may be set to True to not exactly sample the
                        order r of an element g selected uniformly at random
                        from the multiplicative group of the ring of integers
                        modulo N, but rather to sample g and to then
                        heuristically compute r.

      @param verbose  A flag that may be set to True to print intermediary
                      results. Defaults to True.

      @return   This function has no return value. If the test fails, or if some
                other error occurs, an exception is instead raised. """

  # Select factors.
  if verbose:
    print("Sampling N, please wait, this may take a moment...\n");

  factors = [];
  pis = set();

  for i in range(n):
    if e_max == 1:
      ei = 1;
    else:
      ei = int(1 + sample_integer(e_max - 1));

    while True:
      pi = mpz(2 ** (l - 1) + sample_integer(2 ** (l - 1)));

      if not pi.is_prime():
        continue;

      if pi in pis:
        continue;

      break;

    if verbose:
      print("Selected factor " + str(i) + ":", str(pi) + "^" + str(ei));

    factors.append([pi, ei]);

  # Compute N.
  N = prod([pi ** ei for [pi, ei] in factors]);

  if verbose:
    print("\nSelected N = " + str(N) + "\n");

  # Attempt to solve.
  test_solve_j_for_factors_of_N(
    N = N,
    N_factors = factors,
    c_factor = c_factor,
    c_solve = c_solve,
    B_solve = B_solve,
    B_sample = B_sample,
    sample_g = sample_g,
    verbose = verbose);


def test_all(verbose : bool = True) -> None:

  """ @brief  Executes the test suite described in App. A.3 of [E21b].

      [E21b] Ekerå, M.: "On completely factoring any integer efficiently in a
                         single run of an order-finding algorithm".
                        Quantum Inf. Process. 20(6):205 (2021).

      @remark   This convenience function simply iterates over all combinations
                of l in [256, 512, 1024], n in [2, 5, 10, 25], and e_max in
                [1, 2, 3], calling both test_solve_r_for_factors_of_random_N
                and test_solve_j_for_factors_of_random_N() for each combination.

      @param verbose  A flag that may be set to True to print intermediary
                      results. Defaults to True.

      @return   This function has no return value. If the test fails, or if some
                other error occurs, an exception is instead raised. """

  # Setup and start a timer.
  timer = Timer().start();

  for l in [256, 512, 1024]:
    for n in [2, 5, 10, 25]:
      for e_max in [1, 2, 3]:

        if verbose:
          print("\n");

        print("*** Running test for l =", str(l) + ", n =", str(n) +
          ", e_max =", str(e_max) + "...");

        # Solve the chain r -> N.
        test_solve_r_for_factors_of_random_N(
          l = l,
          n = n,
          e_max = e_max,
          verbose = verbose);

        # Solve the whole chain j -> r -> N. This goes slightly beyond [E21b].
        test_solve_j_for_factors_of_random_N(
          l = l,
          n = n,
          e_max = e_max,
          verbose = verbose);

  # Stop the timer.
  timer.stop();

  # All tests have been executed.
  print("\n*** Time required to setup and execute all tests:", timer);
