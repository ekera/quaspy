""" @brief  A module for factoring N completely given the order r of an element
            g selected uniformly at random from the multiplicative group of the
            ring of integers modulo N, where g need not be explicitly
            specified. This by using the algorithm from [E21b].

    This module furthermore contains convenience functions for first solving
    the frequency j yielded by the quantum part of Shor's order-finding
    algorithm [Shor94] for a positive integer multiple r' of r, and for then
    solving r' for the complete factorization of N. This by using the classical
    post-processing algorithms from [E24] in the first step, and the algorithm
    from [E21b] in the second step.

    Finally, this module contains convenience functions for first solving a list
    of n frequencies [j_1, ..., j_n] yielded by n runs of the quantum part of
    Seifert's variation [Seifert01] of Shor's order-finding algorithm [Shor94]
    for a positive integer multiple r' of r, and for then solving r' for the
    complete factorization of N. This by using the classical post-processing
    algorithms from [E24t] (see Sect. 5.4) and [E21], with supporting functions
    from [E24], in the first step, and by using the algorithm from [E21b] in the
    second step.

    Throughout this module, the algorithms are as described in [E24t], [E24],
    [E21] and [E21b]. The notation is also inherited from said works.

    [Shor94] Shor, P.W.: "Algorithms for Quantum Computation: Discrete
                          Logarithms and Factoring".
                         In: Proceedings from FOCS '94, pp. 124–134 (1994).

    [Seifert01] Seifert, J.-P.: "Using fewer qubits in Shor's factorization
                                 algorithm via simultaneous Diophantine
                                 approximation". In: CT-RSA 2001.
                                Springer LNCS 2020, pp. 319–227 (2001).

    [E21] Ekerå, M.: "Quantum algorithms for computing general discrete
                      logarithms and orders with tradeoffs".
                     J. Math. Cryptol. 15(1), pp. 359–407 (2021).

    [E21b] Ekerå, M.: "On completely factoring any integer efficiently in a
                       single run of an order-finding algorithm".
                      Quantum Inf. Process. 20(6):205 (2021).

    [E24] Ekerå, M.: "On the success probability of quantum order finding".
                     ACM Trans. Quantum Comput. 5(2):11 (2024).

    [E24t] Ekerå, M.: "On factoring integers, and computing discrete logarithms
                       and orders, quantumly".
                      PhD thesis, KTH Royal Institute of Technology (2024). """

from enum import Enum;

from gmpy2 import mpz;
from gmpy2 import gcd;
from gmpy2 import powmod;

from .....utils.timer import Timer;
from .....utils.timeout import Timeout;

from .....math.groups import CyclicGroupElement;
from .....math.groups import IntegerModRingMulSubgroupElement;

from .....math.primes import prime_power_product;
from .....math.random import sample_integer;
from .....math.kappa import kappa;

from .....orderfinding.general.postprocessing.ekera import solve_j_for_r;
from .....orderfinding.general.postprocessing.ekera import SolutionMethods;

from .....orderfinding.general.postprocessing.ekera import \
  solve_multiple_j_for_r;
from .....orderfinding.general.postprocessing.ekera import EnumerationOptions;

from .....orderfinding.general.postprocessing.ekera import B_DEFAULT_SOLVE;

from .....math.lattices.lll import LLL_DEFAULT_DELTA;

from .internal.collection import FactorCollection;

class OptProcessCompositeFactors(Enum):

  """ @brief  An enumeration of optimization options for how the solver is
              to select x, and to process composite factors.

      There are three optimization options:

      1. JOINTLY_MOD_N

         Select x uniformly at random from Z_N^*, for N the number to be
         factored, and exponentiate x modulo N to 2^t o.

         This is as described in the algorithm in Sect. 3.2 of [E21b].

      2. JOINTLY_MOD_Np

         Select x uniformly at random from Z_N'^*, for N' the product of all
         pairwise coprime composite factors of N currently stored in the
         collection, and exponentiate x modulo N' to 2^t o.

         This is as above, but with optimizations from Sect. 3.2.1 of [E21b].

      3. SEPARATELY_MOD_Np

         Select x uniformly at random from Z_N'^*, for N' the product of all
         pairwise coprime composite factors of N currently stored in the
         collection. Exponentiate x modulo N' to 2^t o, as N' runs over the
         pairwise coprime composite factors of N currently stored in the
         collection.

         This is as above, with more optimizations from Sect. 3.2.1 of [E21b].

      Note that all three options are equivalent with respect to their ability
      to find non-trivial factors of N. The options differ only in terms of
      their arithmetic complexity. Although several exponentiations may be
      required when the default option is used, this fact is, in general, more
      than compensated for by the fact that the moduli are smaller, leading the
      default option to outperform the other two options.

      For further details, see Sect. 3.2.1 of [E21b].

      [E21b] Ekerå, M.: "On completely factoring any integer efficiently in a
                         single run of an order-finding algorithm".
                        Quantum Inf. Process. 20(6):205 (2021). """

  JOINTLY_MOD_N = 1;

  JOINTLY_MOD_Np = 2;

  SEPARATELY_MOD_Np = 3;


class IncompleteFactorizationException(Exception):

  """ @brief  An exception that is raised to signal an incomplete factorization.

      This occurs only if an iteration or timeout limit has been specified. """

  def __init__(
      self,
      message : str,
      factors : set[int | mpz]):

    """ @brief  Initializes the exception.

        @param message  A message explaining why the exception is raised.
        @param factors  The factors found. """

    super().__init__(message);
    self.factors = factors;


def solve_r_for_factors(
  r : int | mpz,
  N : int | mpz,
  c : int = 1,
  k : int | None = None,
  timeout : int | None | Timeout = None,
  verbose : bool = False,
  opt_split_factors_with_multiplicity : bool = True,
  opt_report_accidental_factors : bool = True,
  opt_abort_early : bool = True,
  opt_square : bool = True,
  opt_exclude_one : bool = True,
  opt_process_composite_factors : OptProcessCompositeFactors =
    OptProcessCompositeFactors.SEPARATELY_MOD_Np) -> set[int]:

  """ @brief  Attempts to factor N completely given the order r of an element g
              selected uniformly at random from the multiplicative group of the
              ring of integers modulo N, where g need not be explicitly
              specified. This by using the algorithm from [E21b].

      Throughout this function, the algorithms are as described in [E21b]
      and [E24]. The notation is also inherited from said works.

      [E21b] Ekerå, M.: "On completely factoring any integer efficiently in a
                         single run of an order-finding algorithm".
                        Quantum Inf. Process. 20(6):205 (2021).

      [E24] Ekerå, M.: "On the success probability of quantum order finding".
                       ACM Trans. Quantum Comput. 5(2):11 (2024).

      @param r   The order r of g, or a positive integer multiple of r.

      @param N   The integer N.

      @param c   A parameter c >= 1 that specifies the maximum size of the
                 missing cm-smooth component in lambda'(N) when solving r for
                 the complete factorization of N. In this context, m is the bit
                 length of N, and cm-smoothness is defined as in [E21b] (and in
                 [E24], see in particular Sect. 4.1.1 in [E24]).

                 As is explained in [E21b], increasing c increases the success
                 probability, at the expense of increasing the runtime.

      @param k   The maximum number of iterations to perform. Defaults to None.

                 If k is set to None, as many iterations as are necessary to
                 completely factor N will be performed. If k is explicitly
                 specified, and the complete factorization of N has not been
                 found after k iterations, an exception of type
                 IncompleteFactorizationException will be raised.

      @param timeout  A timeout after which an IncompleteFactorizationException
                      will be raised and the computation aborted.

                      The timeout may be represented as an integer specifying
                      the timeout in seconds, or as an instance of the Timeout
                      class. May be set to None, as is the default, in which
                      case no timeout is enforced.

      @param verbose  A flag that may be set to True to print intermediary
                      results. Defaults to False.


      All other parameters control various optimizations. They are documented
      below, and in the code. It is recommended to use the defaults.


      @param opt_split_factors_with_multiplicity  A flag that may either be set
                                                  to True (default option) or to
                                                  False.

        When set to True, as is the default, the solver initially tests if
        gcd(r, N) yields a non-trivial factor of N. If so, the non-trivial
        factor is reported.

        When set to False, the aforementioned test is not performed.

        To understand the test, note that if p^e divides N, for e > 1 an integer
        and p a large prime, then p^(e-1) is likely to also divide r. Note
        furthermore that it is relatively inexpensive to split N by computing
        gcd(r, N), compared to exponentiating modulo N.

        For random N void of small factors, it is very unlikely for prime
        factors to occur with multiplicity. For such N, this optimization is
        hence not of any practical use. It is only if N is likely to have
        factors that occur with multiplicity, for some special reason, that
        this optimization is useful in practice. This optimization may also
        yield non-trivial factors if N has small factors. However, such factors
        would typically first be removed, e.g. by trial division, before calling
        upon these more elaborate factoring techniques.

        Note that for N as in the problem instances setup in App. A.3 of [E21b],
        prime factors are intentionally forced to occur with multiplicity with
        high probability (when e_max > 1). This so as to test that such special
        cases are handled correctly by the solver. For such N, this optimization
        is likely to report non-trivial factors. This served as our rationale
        for including it as an option.

        This optimization is described in earlier works in the literature, see
        for instance [GLMS15] for one such description.

        [GLMS15] Grosshans, F., T. Lawson, F. Morain and B. Smith: "Factoring
                 Safe Semiprimes with a Single Quantum Query".
                 ArXiv 1511.04385 (2015).


      @param opt_report_accidental_factors  A flag that may either be set to
                                            True (default option) or to False.

        When set to True, as is the default, the solver reports non-trivial
        factors of N found "by accident" when sampling x (denoted x_j in [E21b])
        uniformly at random from Z_{N'}^*, for N' equal to N, or to the product
        of all pairwise coprime composite factors of N in the factor collection,
        depending on which option is selected for opt_process_composite_factors.
        When set to False, such non-trivial factors are not reported.

        Note that it is very unlikely for non-trivial factors to be found by
        accident if N is void of small prime factors. It is only if small
        factors of N are not first removed, e.g. by trial division, that factors
        are likely to be found by accident. On the other hand, if we do find
        factors by accident, it is only logical to report them. This served as
        our rationale for including this optimization as an option.

        Note furthermore that if N has small factors when this optimization is
        enabled, and if opt_process_composite_factors is furthermore set to
        OptProcessCompositeFactors.JOINTLY_MOD_N, then the same non-trivial
        factor of N may be repeatedly found by accident when sampling. This may
        generate long printouts.

        Also note that factors that are found "by accident" when sampling g
        uniformly at random from Z_N^* are not reported even if
        opt_report_accidental_factors is set to True. This is because such
        factors, if found in practice, would typically affect for which N order
        finding is performed in the first place.


      @param opt_abort_early  A flag that may either be set to True (default
                              option) or to False.

        When set to True, as is the default, the solver computes x^(2^i o) for
        0, 1, ..., min(s, t), for t as in [E21b] and s the least non-negative
        integer such that x^(2^s o) = 1. When set to False, the solver instead
        computes x^(2^i o) for i = 0, 1, ..., t.

        Note that when opt_square (see below) is set to True, the solver first
        computes x^o. It then takes consecutive squares to form x^(2^i o) for
        each i. It follows that the saving incurred by this optimization is
        fairly minor when the opt_square flag is set to True, as it is trivial
        to square one. It is only if the opt_square flag is set to False that
        the saving is substantial.


      @param opt_square   A flag that may either be set to True (default option)
                          or to False.

        When set to True, as is the default, the solver first computes x^o. It
        then takes consecutive squares to form x^(2^i o) for each i. When set to
        False, the solver naïvely computes x^(2^i o) from scratch for each i.

        This optimization is described in Sect. 3.2.1 of [E21b].


      @param opt_exclude_one  A flag that may be set either to True (default
                              option) or to False.

        When set to False, the solver selects g and x (denoted x_j in [E21b])
        uniformly at random from Z_N^* and Z_{N'}^*, respectively, for N' equal
        to N, or to the product of all pairwise coprime composite factors of N
        in the factor collection, depending on which option is selected for the
        opt_process_composite_factors flag.

        When set to True, as is the default, the solver excludes one by instead
        selecting g and x uniformly at random from Z_N^* \ {1} and
        Z_{N'}^* \ {1}, respectively.

        This optimization is described in Sect. 3.2.1 of [E21b].


      @param opt_process_composite_factors  An enumeration entry from the
                                            OptProcessCompositeFactors class
                                            that specifies how x should be
                                            selected, and composite factors
                                            processed, by the solver.

        As is described in OptProcessCompositeFactors, there are three options:

        1. OptProcessCompositeFactors.JOINTLY_MOD_N

          The solver selects x (denoted x_j in [E21b]) uniformly at random from
          Z_N^*. It then exponentiates x modulo N. This is how the unoptimized
          algorithm is described in Sect. 3.2 of [E21b].

        2. OptProcessCompositeFactors.JOINTLY_MOD_Np

          The solver selects x uniformly at random from Z_{N'}^*. It then
          exponentiates x modulo N'.

          Above N' is the product of all pairwise coprime composite factors of N
          currently stored in the factor collection.

        3. OptProcessCompositeFactors.SEPARATELY_MOD_Np (default option)

          The solver selects x uniformly at random from Z_{N'}^*, for N' the
          product of all pairwise coprime composite factors of N currently
          stored in the factor collection.

          The solver then exponentiates x modulo N' where N' now runs over all
          pairwise coprime composite factors of N currently stored in the factor
          collection. This is the default option.

        Note that all three options are equivalent with respect to their ability
        to find non-trivial factors of N. The options differ only in terms of
        their arithmetic complexity. Although several exponentiations may be
        required when the default option is used, this fact is, in general, more
        than compensated for by the fact that the moduli are smaller, leading
        the default option to outperform the other two options.

        This optimization is described in Sect. 3.2.1 of [E21b].


      @return   A set of all distinct prime factors that divide N. """

  # Sanity checks.
  if (r < 1) or (N < 2) or (c < 1):
    raise Exception("Error: Incorrect parameters.");

  # Initial setup.
  timeout = Timeout.parse(timeout);

  # Note: Step 1 is already completed.
  r = mpz(r);
  N = mpz(N);
  m = N.bit_length();

  # Setup and start a timer to measure the total time required to solve.
  timer = Timer().start();

  # Setup and reset a timer to measure the time spent exponentiating.
  timer_exponentiation = Timer();

  # Step 2: Build the product of prime factors q^e < cm and multiply onto r.
  rp = prime_power_product(c * m) * r;

  # Step 3: Let rp = 2^t o for o odd.
  t = kappa(rp);
  o = rp // (2 ** t);

  # Define a pairwise coprime set and add in N.
  F = FactorCollection(N);

  # Optimization: Initially split N when factors of N occur with multiplicity.
  #
  # If p^e divides N for e > 1, then p^(e - 1) is likely to divide r. We may use
  # this fact to initially split N when prime factors occur with multiplicity.
  #
  # Note that splitting N in this way is advantageous, as it can be done without
  # exponentiating, and as it may speed up the subsequent exponentiations (when
  # opt_process_composite_factors is not set to JOINTLY_MOD_N).
  #
  # For further details, see the documentation for solve_r_for_factors().
  if opt_split_factors_with_multiplicity:
    d = gcd(r, N);
    if d != 1:
      if verbose:
        print(
          "Note: Splitting N by gcd(r, N) before commencing to iterate...\n");

      F.add(d);

  # Step 4: For j = 1, 2, ... up to k where k is unbounded.
  j = 0;

  while True:
    # Print current status information before proceeding.
    if verbose:
      print("Iteration:", j);
      F.print_status();

    # Check if we are done...
    if F.is_complete():
      break;

    # Increment j for the next iteration.
    j += 1;

    # Check if j > k, if k is specified, and if so raise an exception passing
    # along the factors that have been found thus far.
    if (k != None) and (j > k):
      raise IncompleteFactorizationException(
        "Error: The iteration limit has been exceeded.", F.found_factors);

    # Check if the timeout has elapsed, if specified, and if so raise an
    # exception passing along the factors that have been found thus far.
    if timeout.is_elapsed():
      raise IncompleteFactorizationException(\
        "Error: The timeout has elapsed.", F.found_factors);

    # Step 4.1: Select x uniformly at random from Z_N^*.

    # Optimization: Select x uniformly at random from Z_N'^*, for N' the product
    # of all pairwise coprime composite factors of N stored in the collection.
    # This as opposed to selecting x uniformly at random from Z_N'^*, for
    # N' = N, when not applying the optimization.
    #
    # For further details, see the documentation for solve_r_for_factors(), and
    # Sect. 3.2.1 of [E21b].
    if opt_process_composite_factors == \
      OptProcessCompositeFactors.JOINTLY_MOD_N:
      # Let N' (denoted Np in the code) be N when not applying the optimization.
      Np = N;
    elif opt_process_composite_factors in \
      [OptProcessCompositeFactors.JOINTLY_MOD_Np,
       OptProcessCompositeFactors.SEPARATELY_MOD_Np]:
      # Let N' (denoted Np in the code) be the product of all pairwise coprime
      # composite factors of N stored in the collection.
      Np = F.residual;
    else:
      raise Exception("Error: Invalid option: opt_process_composite_factors.");

    while True:
      # Sample x uniformly at random from Z_N'^*.
      x = mpz(sample_integer(Np));
      if x == 0:
        continue; # Not in Z_N'^*, and not a non-trivial factor of N'.

      # Optimization: Sample x uniformly at random from Z_N'^* \ {1}.
      #
      # For further details, see the documentation for solve_r_for_factors(),
      # and Sect. 3.2.1 of [E21b].
      if (x == 1) and opt_exclude_one:
        if verbose:
          print("Note: Sampled x = 1; excluding and sampling again...\n");
        continue;

      d = gcd(x, Np);
      if d == 1:
        break; # The element is in Z_N'^*.

      # Optimization: Report the non-trivial factor d found "by accident" above.
      #
      # For further details, see the documentation for solve_r_for_factors().
      if opt_report_accidental_factors:
        if verbose:
          print("Note: Reporting a factor (" + str(d) + ") found \"by " +
                "accident\" when sampling. This is likely to occur only if N " +
                "has small factors.\n");

        # Report the factor.
        F.add(d);

        # Print status again.
        if verbose:
          F.print_status();

        # Check again if we are done, since we have updated F.
        if F.is_complete():
          break;

        if opt_process_composite_factors in \
          [OptProcessCompositeFactors.JOINTLY_MOD_Np,
           OptProcessCompositeFactors.SEPARATELY_MOD_Np]:
          Np = F.residual; # Update to the potentially new N'.

    # Check again if we are done, since we may have updated F above.
    if F.is_complete():
      break;

    # Optimization: Exponentiate x modulo N', for N' the product of all pairwise
    # coprime composite factors of N stored in the collection, or as N' runs
    # over the pairwise coprime composite factor of N stored in the collection.
    # This as opposed to exponentiating x modulo N', for N' = N, when not
    # applying the optimization.
    #
    # For further details, see the documentation for solve_r_for_factors(), and
    # Sect. 3.2.1 of [E21b].
    if opt_process_composite_factors == \
      OptProcessCompositeFactors.SEPARATELY_MOD_Np:
      # Note: For SEPARATELY_MOD_Np, we compute the set of composite pairwise
      # coprime factors stored in the collection and let N' run over the set.
      factors = F.found_factors.difference(F.found_primes);
    elif opt_process_composite_factors in \
      [OptProcessCompositeFactors.JOINTLY_MOD_Np,
       OptProcessCompositeFactors.JOINTLY_MOD_N]:
      # Note: For JOINTLY_MOD_N we have N' = N (where we recall that N' is
      # denoted Np in the code). For JOINTLY_MOD_Np, we have that N' is the
      # product of all pairwise coprime composite factors of N stored in the
      # collection. (Note: This by the manner in which N' was setup above.)
      factors = set([Np]);
    else:
      raise Exception("Error: Invalid option: opt_process_composite_factors.");

    # Exponentiate x for each factor in the set setup above.
    #
    # Note that when opt_process_composite_factors is set to SEPARATELY_MOD_Np,
    # for each composite factor N' processed, any non-trivial factors of N'
    # reported can only split N', as N' is coprime with all other factors of N
    # stored in the collection. This implies that there is no need to go back
    # and re-examine the factor collection after it has been updated with the
    # non-trivial factors reported.
    for Np in factors:
      # Rp = IntegerModRing(Np); # Define the subring Z_N'^*.
      # xp = Rp(x); # Coerce x to Z_N'^*.
      xp = x % Np;

      # Step 4.2: For i = 0, 1, ..., t do:
      timer_exponentiation.start();
      # tmp = xp^o;
      tmp = powmod(xp, o, Np);
      timer_exponentiation.stop();

      # Optimization: Test if tmp = 1, and if so abort early.
      #
      # Were we to proceed, we would obtain d = gcd(tmp^{2^i} - 1, Np) = Np in
      # step 4.2.1, as tmp^{2^i} = 1 for all i, and as gcd(0, Np) = Np.
      #
      # For further details, see the documentation for solve_r_for_factors(),
      # and Sect. 3.2.1 of [E21b].
      if (tmp == 1) and opt_abort_early:
        continue; # No point in continuing with the below procedure.

      # Step 4.2.1 for i = 0.
      d = gcd(tmp - 1, Np);
      if 1 < d < Np:
        F.add(d);

      for i in range(1, t + 1):
        # Optimization: To speed up the arithmetic, we may use a temporary
        # variable tmp, that we initially set to xp^o and then square
        # repeatedly, as opposed to computing xp^(2^i o) for each i.
        #
        # For further details, see the documentation for solve_r_for_factors(),
        # and Sect. 3.2.1 of [E21b].
        timer_exponentiation.start();
        if opt_square:
          tmp = powmod(tmp, 2, Np);
        else:
          tmp = powmod(xp, (2 ** i) * o, Np);
        timer_exponentiation.stop();

        if (tmp == 1) and opt_abort_early:
          break; # No point in continuing to iterate, see above.

        # Step 4.2.1 for i = 1, ..., t.
        d = gcd(tmp - 1, Np);
        if 1 < d < Np:
          F.add(d);

  # Stop the timer.
  timer.stop();

  # The complete factorization has been found.
  if verbose:
    print("Time required to solve:", timer);
    print(" Time spent exponentiating:", timer_exponentiation);
    print(" Time spent checking primality:", F.timer_test_primality);
    print(" Time spent reducing perfect powers:", F.timer_test_perfect_power);

  # Sanity check to assert that the factorization is correct and complete.
  tmp = N;

  for f in F.found_primes:
    if (tmp % f) != 0:
      raise Exception("Error: Failed to factor N correctly.");

    tmp //= f;
    while (tmp % f) == 0:
      tmp //= f;

  if tmp != 1:
    raise Exception("Error: Failed to completely factor N.");

  # Return the set of prime factors found.
  return set([int(q) for q in F.found_primes]);


def solve_j_for_factors(
  j : int | mpz,
  m : int,
  l : int,
  g : CyclicGroupElement,
  N : int | mpz,
  c_solve : int = 1,
  c_factor : int = 1,
  B : int = B_DEFAULT_SOLVE,
  k : int | None = None,
  accept_multiple : bool = True,
  method : SolutionMethods = SolutionMethods.LATTICE_BASED_SHORTEST_VECTOR,
  timeout : int | None | Timeout = None,
  verbose : bool = False,
  opt_speculative : bool = True) -> set[int] | None:

  """ @brief  Attempts to factor N completely, given the frequency j yielded by
              a run of the quantum part of Shor's order-finding algorithm
              [Shor94] when computing the order r of g modulo N, for g an
              element selected uniformly at random from the multiplicative group
              of the ring of integers modulo N.

              This by using the algorithm from [E21b] to factor N given r, or a
              positive multiple of r, and the post-processing algorithm from
              [E24] to find r, or a positive multiple of r, given j.

      Throughout this function, the algorithms are as described in [E24] and
      [E21b]. The notation is also inherited from said works.

      [Shor94] Shor, P.W.: "Algorithms for Quantum Computation: Discrete
                            Logarithms and Factoring".
                           In: Proceedings from FOCS '94, pp. 124–134 (1994).

      [E21b] Ekerå, M.: "On completely factoring any integer efficiently in a
                         single run of an order-finding algorithm".
                        Quantum Inf. Process. 20(6):205 (2021).

      [E24] Ekerå, M.: "On the success probability of quantum order finding".
                       ACM Trans. Quantum Comput. 5(2):11 (2024).

      @remark   This convenience function simply calls solve_j_for_r(), and then
                solve_r_for_factors(), passing along r. To access all options of
                these functions, call them manually in sequence instead.

      @param j  The frequency j yielded by the quantum part of the
                order-finding algorithm.

      @param m  A positive integer m such that r < 2^m.

      @param l  A positive integer l <= m such that m + l is the length of the
                control register in the quantum part of the order-finding
                algorithm.

                If method is set to SolutionMethods.CONTINUED_FRACTIONS_BASED or
                SolutionMethods.LATTICE_BASED_SHORTEST_VECTOR, it is required
                that r^2 < 2^(m + l) or else r may not be found.

                If method is set to SolutionMethods.LATTICE_BASED_ENUMERATE, it
                is possible to select l = m - Delta for some Delta in [0, m), at
                the expense of enumerating at most 6 * sqrt(3) * 2^Delta lattice
                vectors for each offset in j considered.

      @param g  The group element g of order r.

      @param N  The integer N.

      @param c_solve  A parameter c_solve >= 1 that specifies the maximum size
                      of the missing smooth component d in r = d * r_tilde when
                      solving j for r, or a multiple of r. As is explained in
                      [E24], increasing c increases the success probability, at
                      the expense of increasing the runtime.

      @param c_factor   A parameter c_factor >= 1 that specifies the maximum
                        size of the missing smooth component in lambda'(N) when
                        solving r for the complete factorization of N. As is
                        explained in [E21b], increasing c increases the success
                        probability, at the expense of increasing the runtime.

      @param B  A bound B >= 0 on the offset in j. If B > 0, the solve_j_for_r()
                function called by this function tries to solve not only j, but
                also j ± 1, ..., j ± B, for r, or for a positive integer
                multiple of r.

      @param k  The maximum number of iterations to perform when factoring N
                given the order r, or a multiple of r. Defaults to None.

                If k is set to None, as many iterations as are necessary to
                completely factor N will be performed. If k is explicitly
                specified, and the complete factorization of N has not been
                found after k iterations, an exception of type
                IncompleteFactorizationException will be raised.

      @param accept_multiple  A flag that may be set to True to indicate that
                              only a positive integer multiple of r is sought.
                              If set to True, the solve_j_for_r() function
                              called by this function returns as soon as it
                              finds r such that g^r = 1.

      @param method   An enumeration entry from the SolutionMethods class that
                      specifies the method to use to solve j for r. For further
                      details, see the documentation for the SolutionMethods
                      class.

      @param timeout  A timeout after which an IncompleteFactorizationException
                      or TimeoutError will be raised and the computation
                      aborted. More specifically, if the process of factoring N
                      given r, or a multiple of r, has been initiated, an
                      IncompleteFactorizationException will be raised, otherwise
                      a TimeoutError will be raised, when the timeout elapses.

                      The timeout may be represented as an integer specifying
                      the timeout in seconds, or as an instance of the Timeout
                      class. May be set to None, as is the default, in which
                      case no timeout is enforced.

      @param verbose  A flag that may be set to True to print intermediary
                      results and status updates when executing the
                      post-processing algorithm.

      @param opt_speculative  A flag that may be set to True to indicate that
                              Alg. 2 in [E24] should be used instead of Alg. 3
                              to find the missing cm-smooth component of r. In
                              most cases, Alg. 2 is faster than Alg. 3, but in
                              the worst case Alg. 2 is a lot slower than Alg. 3.
                              For further details, see [E24].

      @return   The set of all distinct prime factors that divide N, or None,
                if a positive multiple of r could not be found given j. """

  # Initial setup.
  timeout = Timeout.parse(timeout);
  timeout.check();

  if verbose:
    print("*** Computing a multiple of the order...\n");

  r = solve_j_for_r(
        j = j,
        m = m,
        l = l,
        g = g,
        c = c_solve,
        B = B,
        accept_multiple = accept_multiple,
        method = method,
        timeout = timeout,
        verbose = verbose,
        opt_speculative = opt_speculative);

  if None == r:
    return None;

  if verbose:
    print("\nComputed r = " + str(r) + "\n");

    print("*** Computing the factorization...\n");

  return solve_r_for_factors(
            r = r,
            N = N,
            c = c_factor,
            k = k,
            timeout = timeout,
            verbose = verbose);


def solve_j_for_factors_mod_N(
  j : int | mpz,
  m : int,
  l : int,
  g : int | mpz,
  N : int | mpz,
  c_solve : int = 1,
  c_factor : int = 1,
  B : int = B_DEFAULT_SOLVE,
  k : int | None = None,
  accept_multiple : bool = True,
  method : SolutionMethods = SolutionMethods.LATTICE_BASED_SHORTEST_VECTOR,
  timeout : int | None | Timeout = None,
  verbose : bool = False,
  opt_speculative : bool = True) -> set[int] | None:

  """ @brief  Attempts to factor N completely, given the frequency j yielded by
              a run of the quantum part of Shor's order-finding algorithm
              [Shor94] when computing the order r of g modulo N, for g an
              element selected uniformly at random from the multiplicative group
              of the ring of integers modulo N.

              This by using the algorithm from [E21b] to factor N given r, or a
              positive multiple of r, and the post-processing algorithm from
              [E24] to find r, or a positive multiple of r, given j.

      Throughout this function, the algorithms are as described in [E24] and
      [E21b]. The notation is also inherited from said works.

      [Shor94] Shor, P.W.: "Algorithms for Quantum Computation: Discrete
                            Logarithms and Factoring".
                           In: Proceedings from FOCS '94, pp. 124–134 (1994).

      [E21b] Ekerå, M.: "On completely factoring any integer efficiently in a
                         single run of an order-finding algorithm".
                        Quantum Inf. Process. 20(6):205 (2021).

      [E24] Ekerå, M.: "On the success probability of quantum order finding".
                       ACM Trans. Quantum Comput. 5(2):11 (2024).

      @remark   This convenience function simply calls solve_j_for_factors()
                with g setup by calling IntegerModRingMulSubgroupElement(g, N).

      @remark   In turn, the solve_j_for_factors() convenience function simply
                calls solve_j_for_r(), and then solve_r_for_factors(), passing
                along r. To access all options of these functions, call them
                manually in sequence instead.

      @param j  The frequency j yielded by the quantum part of the
                order-finding algorithm.

      @param m  A positive integer m such that r < 2^m.

      @param l  A positive integer l <= m such that m + l is the length of the
                control register in the quantum part of the order-finding
                algorithm.

                If method is set to SolutionMethods.CONTINUED_FRACTIONS_BASED or
                SolutionMethods.LATTICE_BASED_SHORTEST_VECTOR, it is required
                that r^2 < 2^(m + l) or else r may not be found.

                If method is set to SolutionMethods.LATTICE_BASED_ENUMERATE, it
                is possible to select l = m - Delta for some Delta in [0, m), at
                the expense of enumerating at most 6 * sqrt(3) * 2^Delta lattice
                vectors for each offset in j considered.

      @param g  The group element g of order r.

      @param N  The integer N.

      @param c_solve  A parameter c_solve >= 1 that specifies the maximum size
                      of the missing smooth component d in r = d * r_tilde when
                      solving j for r, or a multiple of r. As is explained in
                      [E24], increasing c increases the success probability, at
                      the expense of increasing the runtime.

      @param c_factor   A parameter c_factor >= 1 that specifies the maximum
                        size of the missing smooth component in lambda'(N) when
                        solving r for the complete factorization of N. As is
                        explained in [E21b], increasing c increases the success
                        probability, at the expense of increasing the runtime.

      @param B  A bound B >= 0 on the offset in j. If B > 0, the
                solve_j_for_r() function called by this function tries
                to solve not only j, but also j ± 1, ..., j ± B, for r,
                or for a positive integer multiple of r.

      @param k  The maximum number of iterations to perform when factoring N
                given the order r, or a multiple of r. Defaults to None.

                If k is set to None, as many iterations as are necessary to
                completely factor N will be performed. If k is explicitly
                specified, and the complete factorization of N has not been
                found after k iterations, an exception of type
                IncompleteFactorizationException will be raised.

      @param accept_multiple  A flag that may be set to True to indicate that
                              only a positive integer multiple of r is sought.
                              If set to True, the solve_j_for_r() function
                              called by this function returns as soon as it
                              finds r such that g^r = 1.

      @param method   An enumeration entry from the SolutionMethods class that
                      specifies the method to use to solve j for r. For further
                      details, see the documentation for the SolutionMethods
                      class.

      @param timeout  A timeout after which an IncompleteFactorizationException
                      or TimeoutError will be raised and the computation
                      aborted. More specifically, if the process of factoring N
                      given r, or a multiple of r, has been initiated, an
                      IncompleteFactorizationException will be raised, otherwise
                      a TimeoutError will be raised, when the timeout elapses.

                      The timeout may be represented as an integer specifying
                      the timeout in seconds, or as an instance of the Timeout
                      class. May be set to None, as is the default, in which
                      case no timeout is enforced.

      @param verbose  A flag that may be set to True to print intermediary
                      results and status updates when executing the
                      post-processing algorithm.

      @param opt_speculative  A flag that may be set to True to indicate that
                              Alg. 2 in [E24] should be used instead of Alg. 3
                              to find the missing cm-smooth component of r. In
                              most cases, Alg. 2 is faster than Alg. 3, but in
                              the worst case Alg. 2 is a lot slower than Alg. 3.
                              For further details, see [E24].

      @return   The set of all distinct prime factors that divide N, or None,
                if a positive multiple of r could not be found given j. """

  g = IntegerModRingMulSubgroupElement(g, N);

  return solve_j_for_factors(
            j = j,
            m = m,
            l = l,
            g = g,
            N = N,
            c_solve = c_solve,
            c_factor = c_factor,
            B = B,
            k = k,
            accept_multiple = accept_multiple,
            method = method,
            timeout = timeout,
            verbose = verbose,
            opt_speculative = opt_speculative);


def solve_multiple_j_for_factors(
  j_list : list[int | mpz],
  m : int,
  l : int,
  g : CyclicGroupElement,
  N : int | mpz,
  tau : int = 0,
  c_solve : int = 1,
  c_factor : int = 1,
  k : int | None = None,
  delta : float = LLL_DEFAULT_DELTA,
  precision : int | None = None,
  enumerate : bool | EnumerationOptions = False,
  timeout : int | None | Timeout = None,
  verbose : bool = False,
  opt_speculative : bool = True) -> set[int] | None:

  """ @brief  Attempts to factor N completely, given a list of n frequencies
              [j_1, ..., j_n] yielded by n independent runs of the quantum part
              of Seifert's variation [Seifert01] of Shor's order-finding
              algorithm [Shor94] when computing the order r of g modulo N, for g
              an element selected uniformly at random from the multiplicative
              group of the ring of integers modulo N.

              This by using the algorithm from [E21b] to factor N given r, or a
              positive multiple of r, and the post-processing algorithm from
              [E24t] (see Sect. 5.4) and [E21], with supporting functions from
              [E24], to find r, or a positive multiple of r, given the n
              frequencies [j_1, ..., j_n].

      Throughout this function, the algorithms are as described in [E24t], [E21]
      and [E21b]. The notation is also inherited from said works.

      [Shor94] Shor, P.W.: "Algorithms for Quantum Computation: Discrete
                            Logarithms and Factoring".
                           In: Proceedings from FOCS '94, pp. 124–134 (1994).

      [Seifert01] Seifert, J.-P.: "Using fewer qubits in Shor's factorization
                                   algorithm via simultaneous Diophantine
                                   approximation". In: CT-RSA 2001.
                                  Springer LNCS 2020, pp. 319–227 (2001).

      [E21] Ekerå, M.: "Quantum algorithms for computing general discrete
                        logarithms and orders with tradeoffs".
                       J. Math. Cryptol. 15(1), pp. 359–407 (2021).

      [E21b] Ekerå, M.: "On completely factoring any integer efficiently in a
                         single run of an order-finding algorithm".
                        Quantum Inf. Process. 20(6):205 (2021).

      [E24] Ekerå, M.: "On the success probability of quantum order finding".
                       ACM Trans. Quantum Comput. 5(2):11 (2024).

      [E24t] Ekerå, M.: "On factoring integers, and computing discrete
                         logarithms and orders, quantumly".
                        PhD thesis, KTH Royal Institute of Technology (2024).

      @remark   This convenience function simply calls solve_multiple_j_for_r(),
                and then solve_r_for_factors(), passing along r. To access all
                options of these functions, call them manually in sequence
                instead.

      @param j_list   The n frequencies [j_1, ..., j_n] where j_1, ..., j_n are
                      integers on [0, 2^(m + l)).

      @param m  A positive integer m such that r < 2^m.

      @param l  A positive integer l ≈ m / s for s a tradeoff factor. The length
                of the control register is m + l in the quantum part of the
                order-finding algorithm.

      @param g  The group element g of order r.

      @param N  The integer N.

      @param tau  An integer tau on [0, l]. Used to scale the basis for the
                  lattice L^tau that is used in the post-processing, and that is
                  generated by the vector (j_1, ..., j_n, 2^tau), and by the n
                  vectors (2^(m + l), 0, ..., 0) thru (0, ..., 0, 2^(m + l), 0).

      @param c_solve  A parameter c_solve >= 1 that specifies the maximum size
                      of the missing smooth component d in r = d * r_tilde when
                      solving j for r, or a multiple of r. As is explained in
                      [E24], increasing c increases the success probability, at
                      the expense of increasing the runtime.

      @param c_factor   A parameter c_factor >= 1 that specifies the maximum
                        size of the missing smooth component in lambda'(N) when
                        solving r for the complete factorization of N. As is
                        explained in [E21b], increasing c increases the success
                        probability, at the expense of increasing the runtime.

      @param k  The maximum number of iterations to perform when factoring N
                given the order r, or a multiple of r. Defaults to None.

                If k is set to None, as many iterations as are necessary to
                completely factor N will be performed. If k is explicitly
                specified, and the complete factorization of N has not been
                found after k iterations, an exception of type
                IncompleteFactorizationException will be raised.

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
                        the lattice L^tau (until r or a positive multiple of r
                        is found or the specified timeout has elapsed), or to
                        EnumerationOptions.SVP to consider only a shortest
                        non-zero vector in L^tau as returned by performing a
                        limited enumeration, or to False to consider only a
                        shortest non-zero vector as returned by LLL.

                        May also be set to EnumerationOptions.BOUNDED_BY_TAU or
                        to EnumerationOptions.BOUNDED_BY_TAU_COMPLETE in which
                        case all vectors within distance R of the origin of the
                        lattice L^tau are enumerated, where R depends on tau as
                        R = sqrt(n + 1) * 2^(m + tau). In the former case, the
                        enumeration is aborted early as soon as r or a positive
                        multiple of r is found. In the latter case, the
                        enumeration runs to completion, after which the minimum
                        candidate for r is taken as r.

      @param timeout  A timeout after which an IncompleteFactorizationException
                      or TimeoutError will be raised and the computation
                      aborted. More specifically, if the process of factoring N
                      given r, or a multiple of r, has been initiated, an
                      IncompleteFactorizationException will be raised, otherwise
                      a TimeoutError will be raised, when the timeout elapses.

                      The timeout may be represented as an integer specifying
                      the timeout in seconds, or as an instance of the Timeout
                      class. May be set to None, as is the default, in which
                      case no timeout is enforced.

      @param verbose  A flag that may be set to True to print intermediary
                      results and status updates when executing the
                      post-processing algorithm.

      @param opt_speculative  A flag that may be set to True to indicate that
                              Alg. 2 in [E24] should be used instead of Alg. 3
                              to find the missing cm-smooth component of r. In
                              most cases, Alg. 2 is faster than Alg. 3, but in
                              the worst case Alg. 2 is a lot slower than Alg. 3.
                              For further details, see [E24].

      @return   The set of all distinct prime factors that divide N, or None,
                if a positive multiple of r could not be found given the n
                frequencies [j_1, ..., j_n]. """

  # Initial setup.
  timeout = Timeout.parse(timeout);
  timeout.check();

  if verbose:
    print("*** Computing a multiple of the order...\n");

  r = solve_multiple_j_for_r(
        j_list =  j_list,
        m = m,
        l = l,
        g = g,
        tau = tau,
        c = c_solve,
        delta = delta,
        precision = precision,
        enumerate = enumerate,
        timeout = timeout,
        verbose = verbose,
        opt_speculative = opt_speculative);

  if None == r:
    return None;

  if verbose:
    print("\nComputed r = " + str(r) + "\n");

    print("*** Computing the factorization...\n");

  return solve_r_for_factors(
            r = r,
            N = N,
            c = c_factor,
            k = k,
            timeout = timeout,
            verbose = verbose);


def solve_multiple_j_for_factors_mod_N(
  j_list : list[int | mpz],
  m : int,
  l : int,
  g : int | mpz,
  N : int | mpz,
  tau : int = 0,
  c_solve : int = 1,
  c_factor : int = 1,
  k : int | None = None,
  delta : float = LLL_DEFAULT_DELTA,
  precision : int | None = None,
  enumerate : bool | EnumerationOptions = False,
  timeout : int | None | Timeout = None,
  verbose : bool = False,
  opt_speculative : bool = True) -> set[int] | None:

  """ @brief  Attempts to factor N completely, given a list of n frequencies
              [j_1, ..., j_n] yielded by n independent runs of the quantum part
              of Seifert's variation [Seifert01] of Shor's order-finding
              algorithm [Shor94] when computing the order r of g modulo N, for g
              an element selected uniformly at random from the multiplicative
              group of the ring of integers modulo N.

              This by using the algorithm from [E21b] to factor N given r, or a
              positive multiple of r, and the post-processing algorithm from
              [E24t] (see Sect. 5.4) and [E21], with supporting functions from
              [E24], to find r, or a positive multiple of r, given the n
              frequencies [j_1, ..., j_n].

      Throughout this function, the algorithms are as described in [E24t], [E21]
      and [E21b]. The notation is also inherited from said works.

      [Shor94] Shor, P.W.: "Algorithms for Quantum Computation: Discrete
                            Logarithms and Factoring".
                           In: Proceedings from FOCS '94, pp. 124–134 (1994).

      [Seifert01] Seifert, J.-P.: "Using fewer qubits in Shor's factorization
                                   algorithm via simultaneous Diophantine
                                   approximation". In: CT-RSA 2001.
                                  Springer LNCS 2020, pp. 319–227 (2001).

      [E21] Ekerå, M.: "Quantum algorithms for computing general discrete
                        logarithms and orders with tradeoffs".
                       J. Math. Cryptol. 15(1), pp. 359–407 (2021).

      [E21b] Ekerå, M.: "On completely factoring any integer efficiently in a
                         single run of an order-finding algorithm".
                        Quantum Inf. Process. 20(6):205 (2021).

      [E24] Ekerå, M.: "On the success probability of quantum order finding".
                       ACM Trans. Quantum Comput. 5(2):11 (2024).

      [E24t] Ekerå, M.: "On factoring integers, and computing discrete
                         logarithms and orders, quantumly".
                        PhD thesis, KTH Royal Institute of Technology (2024).

      @remark   This convenience function simply calls solve_multiple_j_for_r(),
                and then solve_r_for_factors(), passing along r. To access all
                options of these functions, call them manually in sequence
                instead.

      @param j_list   The n frequencies [j_1, ..., j_n] where j_1, ..., j_n are
                      integers on [0, 2^(m + l)).

      @param m  A positive integer m such that r < 2^m.

      @param l  A positive integer l ≈ m / s for s a tradeoff factor. The length
                of the control register is m + l in the quantum part of the
                order-finding algorithm.

      @param g  The group element g of order r.

      @param N  The integer N.

      @param tau  An integer tau on [0, l]. Used to scale the basis for the
                  lattice L^tau that is used in the post-processing, and that is
                  generated by the vector (j_1, ..., j_n, 2^tau), and by the n
                  vectors (2^(m + l), 0, ..., 0) thru (0, ..., 0, 2^(m + l), 0).

      @param c_solve  A parameter c_solve >= 1 that specifies the maximum size
                      of the missing smooth component d in r = d * r_tilde when
                      solving j for r, or a multiple of r. As is explained in
                      [E24], increasing c increases the success probability, at
                      the expense of increasing the runtime.

      @param c_factor   A parameter c_factor >= 1 that specifies the maximum
                        size of the missing smooth component in lambda'(N) when
                        solving r for the complete factorization of N. As is
                        explained in [E21b], increasing c increases the success
                        probability, at the expense of increasing the runtime.

      @param k  The maximum number of iterations to perform when factoring N
                given the order r, or a multiple of r. Defaults to None.

                If k is set to None, as many iterations as are necessary to
                completely factor N will be performed. If k is explicitly
                specified, and the complete factorization of N has not been
                found after k iterations, an exception of type
                IncompleteFactorizationException will be raised.

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
                        the lattice L^tau (until r or a positive multiple of r
                        is found or the specified timeout has elapsed), or to
                        EnumerationOptions.SVP to consider only a shortest
                        non-zero vector in L^tau as returned by performing a
                        limited enumeration, or to False to consider only a
                        shortest non-zero vector as returned by LLL.

                        May also be set to EnumerationOptions.BOUNDED_BY_TAU or
                        to EnumerationOptions.BOUNDED_BY_TAU_COMPLETE in which
                        case all vectors within distance R of the origin of the
                        lattice L^tau are enumerated, where R depends on tau as
                        R = sqrt(n + 1) * 2^(m + tau). In the former case, the
                        enumeration is aborted early as soon as r or a positive
                        multiple of r is found. In the latter case, the
                        enumeration runs to completion, after which the minimum
                        candidate for r is taken as r.

      @param timeout  A timeout after which an IncompleteFactorizationException
                      or TimeoutError will be raised and the computation
                      aborted. More specifically, if the process of factoring N
                      given r, or a multiple of r, has been initiated, an
                      IncompleteFactorizationException will be raised, otherwise
                      a TimeoutError will be raised, when the timeout elapses.

                      The timeout may be represented as an integer specifying
                      the timeout in seconds, or as an instance of the Timeout
                      class. May be set to None, as is the default, in which
                      case no timeout is enforced.

      @param verbose  A flag that may be set to True to print intermediary
                      results and status updates when executing the
                      post-processing algorithm.

      @param opt_speculative  A flag that may be set to True to indicate that
                              Alg. 2 in [E24] should be used instead of Alg. 3
                              to find the missing cm-smooth component of r. In
                              most cases, Alg. 2 is faster than Alg. 3, but in
                              the worst case Alg. 2 is a lot slower than Alg. 3.
                              For further details, see [E24].

      @return   The set of all distinct prime factors that divide N, or None,
                if a positive multiple of r could not be found given the n
                frequencies [j_1, ..., j_n]. """

  g = IntegerModRingMulSubgroupElement(g, N);

  return solve_multiple_j_for_factors(
            j_list = j_list,
            m = m,
            l = l,
            g = g,
            N = N,
            tau = tau,
            c_solve = c_solve,
            c_factor = c_factor,
            k = k,
            delta = delta,
            precision = precision,
            enumerate = enumerate,
            timeout = timeout,
            verbose = verbose,
            opt_speculative = opt_speculative);