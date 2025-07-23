""" @brief  A module for testing the functions in the parent module for
            solving a frequency j yielded by the quantum part of
            Shor's order-finding algorithm for the order r, or optionally for a
            positive integer multiple of r.

    To run all unit tests in this module by calling test_all(), execute:

    ```sh

    $ python3 -m quaspy.orderfinding.general.postprocessing.ekera.test

    ``` """

from ......math.primes import is_B_smooth
from ......utils.timer import Timer;

from ......math.groups import SimulatedCyclicGroupElement;

from ......math.random import sample_integer, sample_l_bit_integer;

from ....sampling import optimal_j_for_z_r, sample_j_given_r;

from ....sampling import B_DEFAULT_SAMPLE;

from ...ekera import B_DEFAULT_SOLVE;

from .. import solve_j_for_r;

from .. import SolutionMethods;

from gmpy2 import mpz;

from math import gcd;

from ..internal.solve import solve_j_for_r_tilde_lattice_enumerate;
from ..internal.solve import solve_j_for_r_tilde_lattice_svp;
from ..internal.solve import solve_j_for_r_tilde_continued_fractions;

def test_solve_j_for_r(
  r : int | mpz,
  m : int | None = None,
  Delta : int | None = None,
  c : int = 1,
  B_solve : int = B_DEFAULT_SOLVE,
  B_sample : int = B_DEFAULT_SAMPLE,
  accept_multiple : bool = False,
  method : SolutionMethods = SolutionMethods.LATTICE_BASED_SHORTEST_VECTOR,
  verbose : bool = True,
  opt_speculative : bool = True,
  opt_isolate_peak : bool = True) -> None:

  """ @brief  Tests the solve_j_for_r() function for a given r, and a given
              set of parameters and optimization flags.

      This function first calls sample_j_given_r(r, m, l, B = B_sample) to
      sample a Fourier frequency j from the distribution induced by the quantum
      order-finding algorithm, for r as passed to this function.

      Unless a larger value of m is explicitly specified, m is set to the bit
      length of the order r. If some Delta in [0, m) is specified, l = m-Delta,
      otherwise l is set to the least positive integer such that r^2 < 2^(m+l).

      @remark   If the sampling fails, this function returns silently, as it is
                expected that sampling j may sometimes fail.

      This function then calls

        solve_j_for_r(
          j,
          m,
          l,
          g,
          c,
          B = B_solve,
          accept_multiple,
          method,
          verbose,
          opt_speculative,
          opt_isolate_peak)

      to attempt to solve j for r, or for some positive multiple of r, depending
      on how the acccept_multiple flag is selected. Above, g is a simulated
      cyclic group element of order r (see the SimulatedCyclicGroupElement
      class), and B = B_solve. All other parameters are simply passed along.

      For the full documentation of these parameters, see solve_j_for_r().

      @param r  The order r.

      @param m  A positive integer m such that r < 2^m, or None, in which case
                this function will set m to the bit length of r.

      @param Delta  An integer in [0, m) used to select l = m-Delta, or None,
                    in which case this function will set l to the least positive
                    integer such that r^2 < 2^(m+l).

      @param c  A parameter c >= 1 that specifies the maximum size of the
                missing cm-smooth component d in r = d * r_tilde when solving j
                for r.

      @param B_solve  A bound B >= 0 on the offset in j when solving j for r.

                      If B > 0, this function tries to solve not only j, but
                      also j ± 1, ..., j ± B, for r, or for a positive integer
                      multiple of r. For further details, see solve_j_for_r().

      @param B_sample   A bound B >= 0 on the offset in j when sampling j. For
                        further details, see sample_j_given_r().


      @param accept_multiple  A flag that may be set to True to indicate that
                              only a positive integer multiple of r is sought.
                              If set to True, this function returns as soon as
                              it finds r such that g^r = 1.

      @param method   An enumeration entry from the SolutionMethods class that
                      specifies the method to use to solve j for r. For further
                      details, see the documentation for the SolutionMethods
                      class.

      @param verbose  A flag that may be set to True to print intermediary
                      results and status updates when executing the test.

      @param opt_speculative  A flag that may be set to True to indicate that
                              Alg. 2 in [E24] should be used instead of Alg. 3
                              to find the missing cm-smooth component of r. In
                              most cases, Alg. 2 is faster than Alg. 3, but in
                              the worst case Alg. 2 is a lot slower than Alg. 3.
                              For further details, see [E24].

      @param opt_isolate_peak   A flag that may be set to True to indicate that
                                all offsets in j up to B should not be tested
                                exhaustively. Rather, the peak around the
                                optimal frequency j_0(z) should be isolated: As
                                soon as offsets to the left and right of j_0(z)
                                are found for which the post-processing
                                algorithm fails to produce r such that g^r = 1,
                                this function returns the minimum r found such
                                that g^r = 1. Note that this flag has no effect
                                if the accept_multiple flag is set to True, as
                                this function then returns as soon as it finds r
                                such that g^r = 1.

      @return   This function has no return value. If the test fails, or if some
                other error occurs, an exception is instead raised.

                The test succeeds if the solve_j_for_r() function behaves as
                expected, i.e. if it succeeds in all cases where it is
                guaranteed to succeed by the analysis in [E24], and either
                succeeds or fails otherwise.

                If the solve_j_for_r() function fails in a case where it is
                guaranteed to succeed in theory, an exception is raised. """

  # Select m and l = m.
  if None == m:
    m = r.bit_length();
  else:
    if r.bit_length() > m:
      raise Exception("Error: It is required that r < 2^m.");

  if None == Delta:
    l = m;
    while (l > 0) and ((r ** 2) < (2 ** (m + l - 1))):
      l -= 1;
  else:
    if (Delta < 0) or (Delta >= m):
      raise Exception("Error: It is required that Delta be on [0, m).");

    l = m - Delta;

  # Sample j.
  sample = sample_j_given_r(r, m, l, B = B_sample, extended_result = True);
  if None == sample:
    if verbose:
      print("\n*** Failed to sample j.\n");

    return;
  else:
    [j, [z, j0, j0_offset]] = sample;

  if verbose:
    print("Sampled j = " + str(j) + "\n");

  # Setup g.
  g = SimulatedCyclicGroupElement(r);

  # Setup and start a timer.
  timer = Timer().start();

  # Attempt to recover r.
  candidate_r = \
    solve_j_for_r(
      j = j,
      m = m,
      l = l,
      g = g,
      c = c,
      B = B_solve,
      accept_multiple = accept_multiple,
      method = method,
      verbose = verbose,
      opt_speculative = opt_speculative,
      opt_isolate_peak = opt_isolate_peak);

  # Stop the timer.
  timer.stop();

  if None == candidate_r:
    # Check if we expect to fail since the offset is too large.
    if abs(j0_offset) > B_solve:
      if verbose:
        print("\n*** Failed to recover r, as expected, " +
          "since the offset in j is greater than B.\n");

      return;

    # Check if we expect to fail since d = gcd(r, z) is not cm-smooth.
    d = gcd(z, r);

    if not is_B_smooth(d, c * m):
      if verbose:
        print("\n*** Failed to recover r, as expected, " +
          "since d = gcd(r, z) is not cm-smooth.\n");

      return;

    # Unexpected failure.
    if verbose:
      print("\n*** Failed to recover r.\n");

    raise Exception("Error: Failed to solve for r, despite the problem " +
      "instance being solvable in theory.");

  # Compare.
  if candidate_r == r:
    if verbose:
      print("\n*** Successfully recovered the order r.");
      print("\nTotal time required to solve j for r: " + str(timer) + "\n");

    return;

  if accept_multiple and (candidate_r >= 1) and ((candidate_r % r) == 0):
    if verbose:
      print("\n*** Successfully recovered a multiple of the order r.");
      print("\nTotal time required to solve j for a multiple r: " + \
        str(timer) + "\n");

    return;

  if verbose:
    print("\n*** Failed to recover r.\n");

  raise Exception("Error: Failed to solve for r, despite the problem " +
    "instance being solvable in theory.");


def test_solve_j_for_random_r(
  m : int = 2048,
  Delta : int | None = None,
  B_solve : int = B_DEFAULT_SOLVE,
  B_sample : int = B_DEFAULT_SAMPLE,
  accept_multiple : bool = False,
  method : SolutionMethods = SolutionMethods.LATTICE_BASED_SHORTEST_VECTOR,
  verbose : bool = True,
  opt_speculative : bool = True,
  opt_isolate_peak : bool = True) -> None:

  """ @brief  Selects an m bit order r uniformly at random from the set of all
              such orders, and then calls test_solve_j_for_r() for this order r,
              passing along the other arguments.

      @remark   This is a convenience function. For details on all other
                parameters, please see the documentation of the
                test_solve_j_for_r() function.

      @return   This function has no return value. If the test fails, or if some
                other error occurs, an exception is instead raised. """

  # Select an m-bit order r.
  r = sample_l_bit_integer(m);

  if verbose:
    print("Sampled r = " + str(r) + "\n");

  # Attempt to solve.
  test_solve_j_for_r(
    r = r,
    m = m,
    Delta = Delta,
    B_solve = B_solve,
    B_sample = B_sample,
    accept_multiple = accept_multiple,
    method = method,
    verbose = verbose,
    opt_speculative = opt_speculative,
    opt_isolate_peak = opt_isolate_peak);


def test_all_solve_for_r(verbose : bool = True) -> None:

  """ @brief  Calls test_solve_j_for_random_r() for all combinations of
              accept_multiple and opt_speculative in {True, False}, and, for
              each such combination, for all supported solution methods, and for
              each such combination, for all values of m in
              {128, 256, 384, 512, 1024, 2048, 4096, 8192}, passing along the
              verbose flag.

      More specifically, this function calls test_solve_j_for_random_r() for all
      combinations of accept_multiple and opt_speculative in {True, False}, and,
      for each such combination, for all solution methods in

        {SolutionMethods.CONTINUED_FRACTIONS_BASED,
         SolutionMethods.LATTICE_BASED_SHORTEST_VECTOR;
         SolutionMethods.LATTICE_BASED_ENUMERATE},

      and, for each such combination, for all values of m in

        {128, 256, 384, 512, 1024, 2048, 4096, 8192}.

      If test_solve_j_for_random_r() raises an exception due to a test failing,
      this exception is not caught but simply allowed to propagate.

      This function always sets opt_isolate_peak to True, as is the default.

      @param verbose  A flag that may be set to True to print intermediary
                      results. Defaults to False.

      @return   This function has no return value. If the test fails, or if some
                other error occurs, an exception is instead raised. """

  for accept_multiple in [True, False]:
    for method in [SolutionMethods.CONTINUED_FRACTIONS_BASED,
                   SolutionMethods.LATTICE_BASED_SHORTEST_VECTOR,
                   SolutionMethods.LATTICE_BASED_ENUMERATE]:

      if not verbose:
        print("");

      for opt_speculative in [True, False]:
        for opt_isolate_peak in [True]:

          for m in [128, 256, 384, 512, 1024, 2048, 4096, 8192]:

            if verbose:
              print("");

            if method == SolutionMethods.LATTICE_BASED_ENUMERATE:
              Deltas = [0, 1, 2, 3, 4, 5, 8, 12, 16];
            else:
              Deltas = [None];

            for Delta in Deltas:
              print("*** Running test for " +
                "m =", str(m) + ", " +
                "Delta =", str(Delta) + ", " +
                "accept_multiple =", str(accept_multiple) + ", " +
                "method =", str(method) + ", " +
                "opt_speculative =", str(opt_speculative) + ", " +
                "opt_isolate_peak =", str(opt_isolate_peak) + "...");

              test_solve_j_for_random_r(
                m = m,
                Delta = Delta,
                accept_multiple = accept_multiple,
                method = method,
                verbose = verbose,
                opt_speculative = opt_speculative,
                opt_isolate_peak = opt_isolate_peak);


def test_solve_for_r_tilde(
  z : int | mpz,
  r : int | mpz,
  m : int,
  c : int = 1,
  Deltas : bool = False,
  verbose : bool = False) -> None:

  """ @brief  Calls solve_j_for_r_tilde_continued_fractions(),
              solve_j_for_r_tilde_lattice_svp() and
              solve_j_for_r_tilde_lattice_enumerate() with j = j0(z) to test
              that r_tilde = r / gcd(r, z) is successfully recovered.

      This function picks the least positive integer l such that r^2 < 2^(m+l),
      computes the optimal frequency j = j0(z) given r, m and l, and attempts
      to solve j for r using

        - the continued fractions-based solver, by calling the
          solve_j_for_r_tilde_continued_fractions() function,

        - the lattice-based non-zero shortest vector problem (SVP) solver, by
          calling the solve_j_for_r_tilde_lattice_svp() function, and

        - the lattice-based solver with enumeration, by calling the
          solve_j_for_r_tilde_lattice_enumerate() function.

      For each solver, it is verified that r_tilde = r / gcd(r, z) is
      successfully recovered.

      If the Deltas flag is set to true, this function also attempt to solve j
      for r using the lattice-based solver with enumeration, whilst setting
      l = m - Delta och re-computing j = j0(z) given r, m and l, for all Delta
      in [0, ..., min(m - 1, 12)]. For each value of Delta, it is again verified
      that r_tilde is successfully  recovered.

      @remark   This function checks whether d = gcd(r, z) is cm-smooth.

                If d is not cm-smooth, all tests of the lattice-based solver
                with support for enumerating the lattice are skipped over the
                solver in question is not guaranteed to find r_tilde unless d is
                cm-smooth.

      @param z  An integer z on [0, r) identifying the peak index.

      @param r  The order r.

      @param m  A positive integer m such that r < 2^m.

      @param c  A parameter c >= 1 that specifies the maximum size of the
                missing cm-smooth component d in r = d * r_tilde.

      @param Deltas   A flag that may be set to True to perform extended tests
                      for the lattice-based solver with enumeration, as
                      explained above. Defaults to False.

      @param verbose  A flag that may be set to True to print intermediary
                      results. Defaults to False.

      @return   This function has no return value. If the test fails, or if some
                other error occurs, an exception is instead raised."""

  l = m;
  while (l > 1) and ((r ** 2) < (2 ** (m + l - 1))):
    l -= 1;

  j0 = optimal_j_for_z_r(z = z, r = r, m = m, l = l);

  d = gcd(r, z);

  smooth = is_B_smooth(d, c * m);

  r_tilde = r // d;

  # Continued fractions:
  r_tilde_candidates = solve_j_for_r_tilde_continued_fractions(j0, m, l);

  if r_tilde != r_tilde_candidates:
    raise Exception(
      "Error: Failed to recover r_tilde by continued fractions.");

  # Lattice SVP:
  [r_tilde_candidates, _] = solve_j_for_r_tilde_lattice_svp(j0, m, l);

  if r_tilde != r_tilde_candidates:
    raise Exception(
      "Error: Failed to recover r_tilde by solving the SVP in L.");

  # Lattice enumeration:
  if smooth:
    g = SimulatedCyclicGroupElement(r);

    [r_tilde_candidates, _] = solve_j_for_r_tilde_lattice_enumerate(
                                j = j0,
                                m = m,
                                l = l,
                                g = g,
                                c = c,
                                accept_multiple = False,
                                filtered_r_tilde_candidates = set(),
                                verbose = verbose);

    if type(r_tilde_candidates) != type(set()):
      raise Exception("Error: Expected a set to the returned.");

    if not (r_tilde in r_tilde_candidates):
      raise Exception("Error: Failed to recover r_tilde by enumerating L.");


    [r_tilde_candidates, _] = solve_j_for_r_tilde_lattice_enumerate(
                                j = j0,
                                m = m,
                                l = l,
                                g = g,
                                c = c,
                                accept_multiple = False,
                                verbose = verbose);

    if not (r_tilde in r_tilde_candidates):
      raise Exception("Error: Failed to recover r_tilde by enumerating L.");


    if Deltas:

      for Delta in range(min(m - 1, 12) + 1):
        l = m - Delta;

        j0 = optimal_j_for_z_r(z = z, r = r, m = m, l = l);


        [r_tilde_candidates, _] = solve_j_for_r_tilde_lattice_enumerate(
                                    j = j0,
                                    m = m,
                                    l = l,
                                    g = g,
                                    c = c,
                                    accept_multiple = False,
                                    filtered_r_tilde_candidates = set(),
                                    verbose = verbose);

        if type(r_tilde_candidates) != type(set()):
          raise Exception("Error: Expected a set to the returned.");

        if not (r_tilde in r_tilde_candidates):
          raise Exception("Error: Failed to recover r_tilde by enumerating L.");


        [r_tilde_candidates, _] = solve_j_for_r_tilde_lattice_enumerate(
                                    j = j0,
                                    m = m,
                                    l = l,
                                    g = g,
                                    c = c,
                                    accept_multiple = False,
                                    verbose = verbose);

        if not (r_tilde in r_tilde_candidates):
          raise Exception("Error: Failed to recover r_tilde by enumerating L.");


def test_solve_for_r_tilde_exhaustive(
  m : int = 16,
  c : int = 1,
  verbose : bool = False) -> None:

  """ @brief  Selects r uniformly at random from the set of all m-bit integers,
              and then calls test_solve_for_r_tilde() for r, and for all values
              of z in [0, r), with Deltas = False, passing along the verbose
              flag.

      @param m  A positive integer m.

      @param c  A parameter c >= 1 that specifies the maximum size of the
                missing cm-smooth component d in r = d * r_tilde.

      @param verbose  A flag that may be set to True to print intermediary
                      results. Defaults to False.

      @return   This function has no return value. If the test fails, or if some
                other error occurs, an exception is instead raised. """

  # Pick an m-bit order r at random.
  r = sample_l_bit_integer(m);

  print("*** Running r_tilde recovery test for r = " + str(r) + "...\n");

  g = SimulatedCyclicGroupElement(r);

  for z in range(r):
    if verbose:
      print(" >> Processing " + str(z) + " / " + str(r));

    test_solve_for_r_tilde(
      z = z,
      r = r,
      m = m,
      c = c,
      Deltas = False,
      verbose = verbose);


def test_solve_for_r_tilde_randomized(
  m : int = 2048,
  c : int = 1,
  verbose : bool = False) -> None:

  """ @brief  Selects r uniformly at random from the set of all m-bit integers,
              then selects z uniformly at random from [0, r), and finally calls
              test_solve_for_r_tilde() for r and z, with Deltas = True, passing
              along the other arguments.

      @param m  A positive integer m.

      @param c  A parameter c >= 1 that specifies the maximum size of the
                missing cm-smooth component d in r = d * r_tilde.

      @param verbose  A flag that may be set to True to print intermediary
                      results. Defaults to False.

      @return   This function has no return value. If the test fails, or if some
                other error occurs, an exception is instead raised. """

  # Pick an m-bit order r at random.
  r = sample_l_bit_integer(m);

  print("*** Running r_tilde recovery test for r = " + str(r) + "...\n");

  z = sample_integer(r);

  test_solve_for_r_tilde(
    z = z,
    r = r,
    m = m,
    c = c,
    Deltas = True,
    verbose = verbose);


def test_all_solve_for_r_tilde(
  verbose : bool = False) -> None:

  """ @brief  Calls test_solve_for_r_tilde_exhaustive(m, verbose) for all
              integers m in {4, 5, ..., 16}, and then calls
              test_solve_for_r_tilde_randomized(m, verbose) ten times
              for m in {256, 512, 1024, 2048, 4096, 8192, 16384}, passing
              along the verbose flag.

      @param verbose  A flag that may be set to True to print intermediary
                      results. Defaults to False.

      @return   This function has no return value. If the test fails, or if some
                other error occurs, an exception is instead raised. """

  for m in range(4, 16 + 1):
    test_solve_for_r_tilde_exhaustive(m = m, verbose = verbose);

  for m in [256, 512, 1024, 2048, 4096, 8192, 16384]:
    for _ in range(10):
      test_solve_for_r_tilde_randomized(m = m, verbose = verbose);


def test_all(
  verbose : bool = False) -> None:

  """ @brief  Calls test_all_solve_for_r_tilde(), and then
              test_all_solve_for_r(), passing along the verbose flag.

      @param verbose  A flag that may be set to True to print intermediary
                      results. Defaults to False.

      @return   This function has no return value. If the test fails, or if some
                other error occurs, an exception is instead raised. """

  # Setup and start a timer.
  timer = Timer().start();

  test_all_solve_for_r_tilde(verbose = verbose);

  test_all_solve_for_r(verbose = verbose);

  # Stop the timer.
  timer.stop();

  if not verbose:
    print("");

  print("*** Time required to setup and execute all tests:", timer);
