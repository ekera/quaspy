""" @brief  A module for solving a frequency j yielded by the quantum part of
            Shor's order-finding algorithm for the order r, or optionally for a
            positive integer multiple of r. This by using the post-processing
            algorithms in [E22p].

    [E22p] Ekerå, M.: "On the success probability of quantum order finding".
                      ArXiv 2201.07791 (2022). """

from enum import Enum;

from gmpy2 import mpz;

from gmpy2 import gcd;

from .....math.primes import prime_power_product;

from .....math.groups import CyclicGroupElement;
from .....math.groups import IntegerModRingMulSubgroupElement;

from .internal.solve import solve_j_for_r_tilde_continued_fractions;
from .internal.solve import solve_j_for_r_tilde_lattice_svp;
from .internal.solve import solve_j_for_r_tilde_lattice_enumerate;

from .internal.collection import CandidateCollection;

from .internal.algorithms import algorithm2;
from .internal.algorithms import algorithm3;

from .. import B_DEFAULT_SOLVE;

class SolutionMethods(Enum):

  """ @brief  An enumeration of methods for solving j for r.

      As is explained in [E22p], these methods are all designed to solve an
      optimal frequency j = j0(z), where z in [0, r), for r_tilde = r / d,
      where d = gcd(r, z) and d is cm-smooth. In what follows below, it is
      assumed that these requirements are met.

      There are three solution methods:

      1. CONTINUED_FRACTIONS_BASED

         Expand j / 2^(m+l) in continued fractions to find z / r, and hence
         r_tilde = r / gcd(r, z), as originally proposed in [Shor94], but with
         slightly smaller m+l, as decribed in [E22p].

         By Lemma 6 in [E22p], the last convergent p / q with denominator
         q < 2^((m+l)/2) in the continued fractions expansion of j / 2^(m+l) is
         equal to z / r, provided that r < 2^m and r^2 < 2^(m+l).

         For further details, see Lemma 6, and Sect. 4 and App. B, of [E22p].

      2. LATTICE_BASED_SHORTEST_VECTOR

         Use Lagrange's lattice basis reduction algorithm to find the shortest
         non-zero vector in the lattice L spanned by (j, 1/2) and (2^(m+l), 0).

         By Lemma 7 in [E22p], provided that r < 2^m and r^2 < 2^(m+l), the
         second component of the shortest non-zero vector in L has r_tilde / 2
         in its second component, up to sign of course.

         For further details, see Lemma 7, and Sect. 4 and App. C, of [E22p].

      3. LATTICE_BASED_ENUMERATE

         Use Lagrange's lattice basis reduction algorithm to find a reduced
         basis for the lattice L spanned by (j, 1/2) and (2^(m+l), 0), and
         enumerate all vectors within a ball of radius 2^(m-1/2) in L centered
         at the origin to find u = (rj - 2^(m+l) z, r / 2) / d, and hence
         r_tilde, as the second component is r_tilde / 2.

         By Lemma 8 in [E22p], provided that r < 2^m and l = m - Delta, at most
         6 * sqrt(3) * 2^Delta vectors must be enumerated in L to find u and
         hence r_tilde, so if Delta is small then this method is efficient.

         In practice, as mentioned in [E22p], the leading constant in the above
         bound is not tight, and the enumeration can be optimized. Some of these
         optimizations are implemented here so the enumeration typically
         considers fewer vectors than the bound indicates.

         For further details, see Lemma 8, and Sect. 4 and App. C, of [E22p].

      [Shor94] Shor, P.W.: "Algorithms for Quantum Computation: Discrete
                            Logarithms and Factoring".
                           In: Proceedings from FOCS '94, pp. 124–134 (1994).

      [E22p] Ekerå, M.: "On the success probability of quantum order finding".
                        ArXiv 2201.07791 (2022). """

  CONTINUED_FRACTIONS_BASED = 1;

  LATTICE_BASED_SHORTEST_VECTOR = 2;

  LATTICE_BASED_ENUMERATE = 3;


def solve_j_for_r(
  j : int,
  m : int,
  l : int,
  g: CyclicGroupElement,
  c : int = 1,
  B : int = B_DEFAULT_SOLVE,
  accept_multiple = False,
  method = SolutionMethods.LATTICE_BASED_SHORTEST_VECTOR,
  verbose = False,
  opt_speculative = True,
  opt_isolate_peak = True):

  """ @brief  Attempts to compute the order r of g, or a positive integer
              multiple thereof, given a frequency j yielded by the quantum part
              of Shor's order-finding algorithm, by using the post-processing
              algorithms described in detail in [E22p].

      [E22p] Ekerå, M.: "On the success probability of quantum order finding".
                        ArXiv 2201.07791 (2022).

      The idea is to try to solve not only j, but also j ± 1, .., j ± B, for r,
      with the aim of solving an optimal frequency j0(z) for r, for z the peak
      index on [0, r). Provided

        - that j0(z) is solved for r,

        - that d = gcd(r, z) is cm-smooth
          (for the definition of cm-smooth in [E22p]),

        - that l is selected as required by the solution method
          (see the documentation for parameter l below), and

        - that the accept_multiple flag is set to False,

      the order r will be found by this function.

      Note that this function does not implement meet-in-the-middle-techniques,
      although it is noted in [E23p] that it is possible to use such techniques
      to speed up the post-processing also for order finding.

      [E23p] Ekerå, M.: "On the success probability of the quantum algorithm for
                         the short DLP". ArXiv 2309.01754 (2023).

      @param j  The frequency j yielded by the quantum order-finding algorithm.

      @param m  A positive integer m such that r < 2^m.

      @param l  A positive integer l <= m, such that m+l is the length of the
                control register in the quantum order-finding algorithm.

                If method is set to SolutionMethods.CONTINUED_FRACTIONS_BASED or
                SolutionMethods.LATTICE_BASED_SHORTEST_VECTOR, it is required
                that r^2 < 2^(m+l) or else r may not be found.

                If method is set to SolutionMethods.LATTICE_BASED_ENUMERATE, it
                is possible to select l = m - Delta for some Delta in [0, m), at
                the expense of enumerating at most 6 * sqrt(3) * 2^Delta lattice
                vectors for each offset in j considered.

      @param g  The group element g of order r.

      @param c  A parameter c >= 1 that specifies the maximum size of the
                missing cm-smooth component d in r = d * r_tilde when solving j
                for r, for cm-smooth as defined in [E22p].

                As is explained in [E22p], increasing c increases the success
                probability, at the expense of increasing the runtime.

      @param B  A bound B >= 0 on the offset in j. If B > 0, this function tries
                to solve not only j, but also j ± 1, .., j ± B, for r, or for a
                positive integer multiple of r.

      @param accept_multiple  A flag that may be set to True to indicate that
                              only a positive integer multiple of r is sought.
                              If set to True, this function returns as soon as
                              it finds r such that g^r = 1.

      @param method   An enumeration entry from the SolutionMethods class that
                      specifies the method to use to solve j for r. For further
                      details, see the documentation for the SolutionMethods
                      class.

      @param verbose  A flag that may be set to True to print intermediary
                      results and status updates when executing the
                      post-processing algorithm.

      @param opt_speculative  A flag that may be set to True to indicate that
                              Algorithm 2 in [E22p] should be used instead of
                              Algorithm 3 to find the missing cm-smooth
                              component of r. In most cases, Algorithm 2 is
                              faster than Algorithm 3, but in the worst case
                              Algorithm 2 is a lot slower than Algorithm 3. For
                              further details, see [E22p].

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

      @return   If the accept_multiple flag is set to False, the order r of g
                is returned with probability >= P, for P as given by the lower
                bound in [E22p] (provided that the opt_isolate_peak flag is set
                to False). Otherwise, None, or exceptionally a positive integer
                multiple of the order r, is returned. If the accept_multiple
                flag is set to True, some positive integer multiple of r is
                returned with probabilty >= P. Otherwise, None is returned. """

  # Configure algorithms.
  if opt_speculative:
    algorithm = algorithm2;
  else:
    algorithm = algorithm3;

  # Setup the set S', and keep it reduced by using a special class.
  filtered_r_tilde_candidates = CandidateCollection();

  # Setup a set of dismissed reduced candidates for r_tilde.
  dismissed_reduced_r_tilde_candidates = set();

  # Initially set mu to zero. For each r_tilde_candidate to be tested, check if
  #
  #   gp^(reduced_r_tilde_candidate) = g^(reduced_r_tilde_candidate * e) = 1
  #
  # for reduced_r_tilde_candidate = gcd(r_tilde_candidate, mu). If so, we know
  # that reduced_r_tilde_candidate * e is a multiple of r, so we can update mu
  # by mu <- gcd(mu, reduced_r_tilde_candidate * e). This also explains why for
  # each r_tilde_candidate, it is sufficient to test reduced_r_tilde_candidate.

  mu = mpz(0);

  # Pre-compute e and gp.
  e = prime_power_product(c * m);

  gp = g ** e;
  if gp == 1:
    # Trivial case.
    return algorithm(g, 1, m, c);

  j = mpz(j);

  pow2ml = 2 ** (m + l);
  if B > pow2ml // 2:
    B = pow2ml // 2;
  pow2ml = mpz(pow2ml);

  [left, right] = [None, None];
  [skip_left, skip_right] = [False, False];

  positive_multiples = None;
  negative_multiples = None;

  for offset in range(0, B + 1):
    if verbose:
      print("Trying offset:", offset, "/", B);

    if skip_left and skip_right:
      break;

    for sign in [1, -1]:
      if (0 == offset) and (-1 == sign):
        continue;

      if opt_isolate_peak:
        if [left, right] != [None, None]:
          if sign == 1:
            if offset > right + 1:
              skip_right = True;
              continue;
          else:
            if -offset < left - 1:
              skip_left = True;
              continue;

      offset_j = (j + sign * offset) % pow2ml;

      # Use the prescribed solution method to find candidates for r_tilde.
      if method == SolutionMethods.LATTICE_BASED_SHORTEST_VECTOR:

        # Fetch an initial guess for the row multiples.
        if offset == 0:
          multiples = None;
        else:
          if sign == 1:
            multiples = positive_multiples;
          else:
            multiples = negative_multiples;

        # Assume that we are not able to solve.
        success = False;

        # Solve for r_tilde.
        [r_tilde_candidates, multiples] = \
          solve_j_for_r_tilde_lattice_svp(offset_j, m, l, multiples);
        r_tilde_candidates = [r_tilde_candidates];

        # Update the guess for the row multiples.
        if offset == 0:
          positive_multiples = negative_multiples = multiples;
        else:
          if sign == 1:
            positive_multiples = multiples;
          else:
            negative_multiples = multiples;

      elif method == SolutionMethods.LATTICE_BASED_ENUMERATE:

        # Fetch an initial guess for the row multiples.
        if offset == 0:
          multiples = None;
        else:
          if sign == 1:
            multiples = positive_multiples;
          else:
            multiples = negative_multiples;

        # Solve for r_tilde.
        [filtered_r_tilde_candidates,
          [success,
           dismissed_reduced_r_tilde_candidates,
           mu,
           multiples]] = \
          solve_j_for_r_tilde_lattice_enumerate(
            j = offset_j,
            m = m,
            l = l,
            g = g,
            g_pow_e_context = [gp, e],
            c = c,
            accept_multiple = accept_multiple,
            filtered_r_tilde_candidates = filtered_r_tilde_candidates,
            dismissed_reduced_r_tilde_candidates =
              dismissed_reduced_r_tilde_candidates,
            multiples = multiples,
            verbose = verbose);

        # Update the guess for the row multiples.
        if offset == 0:
          positive_multiples = negative_multiples = multiples;
        else:
          if sign == 1:
            positive_multiples = multiples;
          else:
            negative_multiples = multiples;

      elif method == SolutionMethods.CONTINUED_FRACTIONS_BASED:

        # Assume that we are not able to solve.
        success = False;

        # Solve for r_tilde.
        r_tilde_candidates = \
          [solve_j_for_r_tilde_continued_fractions(offset_j, m, l)];

      else:
        raise Exception("Error: Incorrect parameters: Unknown method.");


      if method in [SolutionMethods.LATTICE_BASED_SHORTEST_VECTOR,
                    SolutionMethods.CONTINUED_FRACTIONS_BASED]:

        # When using either of these solution methods, the candidates in
        # r_tilde_candidates need not fulfill the requirement that
        #
        #   gp^r_tilde_candidate == 1,
        #
        # so we need to check this requirement. We do this by first checking
        # if r_tilde_candidate is already in the filtered set of candidates for
        # r_tilde that have passed this requirement. Only if r_tilde_candidate
        # is not in filtered_r_tilde_candidates do we need to exponentiate.

        for r_tilde_candidate in r_tilde_candidates:
          if r_tilde_candidate in filtered_r_tilde_candidates:
            # We already found this candidate, or a divisor of it.
            success = True;
          else:
            # Use that mu is an r-multiple to reduce the candidate for r_tilde.
            reduced_r_tilde_candidate = gcd(r_tilde_candidate, mu);

            # Check if the remainder after reduction is equal to one, or if we
            # already have dismissed this reduced candidate. Otherwise proceed.
            if (reduced_r_tilde_candidate == 1) or \
               (reduced_r_tilde_candidate in \
                dismissed_reduced_r_tilde_candidates):
              # Dismiss the reduced candidate.
              if verbose:
                print("Dismissing:", r_tilde_candidate);
            else:
              # The reduced candidate has not already been dismissed.
              if verbose:
                print("Testing the reduced candidate:", \
                  reduced_r_tilde_candidate);

              # Test the reduce candidate.
              if (gp ** reduced_r_tilde_candidate) == 1:
                if accept_multiple:
                  # Return immediately.
                  return algorithm(g, r_tilde_candidate, m, c);

                # Add r_tilde_candidate to the filtered candidates for r_tilde.
                filtered_r_tilde_candidates.add(r_tilde_candidate);

                # We know that reduced_r_tilde_candidate * e is a multiple of r,
                # so we may update mu to account for this fact:
                mu = gcd(reduced_r_tilde_candidate * e, mu);

                # Mark this j as successful.
                success = True;
              else:
                # Add reduced_r_tilde_candidate to the dismissed reduced
                # candidates for r_tilde to avoid future exponentiations.
                dismissed_reduced_r_tilde_candidates.\
                  add(reduced_r_tilde_candidate);

      elif method == SolutionMethods.LATTICE_BASED_ENUMERATE:

        # When using this solution method, the enumerate_lattice() function will
        # already have tested the requirement that
        #
        #   gp^r_tilde_candidate == 1,
        #
        # as very many candidates could otherwise be returned. Hence, there is
        # no reason to repeat the test here by exponentiating.

        if accept_multiple:
          if len(filtered_r_tilde_candidates) != 0:
            # Return immediately.
            return algorithm(g, min(filtered_r_tilde_candidates), m, c);

      else:
        raise Exception("Error: Incorrect parameters: Unknown method.");

      if opt_isolate_peak and success:
        if [left, right] == [None, None]:
          left = right = offset;
        elif sign == 1:
          if right == offset - 1:
            right = offset;
          else:
            raise Exception("Error: Not a contiguous interval for the peak.");
        else:
          if left == -offset + 1:
            left = -offset;
          else:
            raise Exception("Error: Not a contiguous interval for the peak.");

  if len(filtered_r_tilde_candidates) == 0:
    return None;

  return min([algorithm(g, r_tilde_candidate, m, c)
    for r_tilde_candidate in filtered_r_tilde_candidates])


def solve_j_for_r_mod_N(
  j : int,
  m : int,
  l : int,
  g: int,
  N: int,
  c : int = 1,
  B : int = B_DEFAULT_SOLVE,
  accept_multiple = False,
  method = SolutionMethods.LATTICE_BASED_SHORTEST_VECTOR,
  verbose = False,
  opt_speculative = True,
  opt_isolate_peak = True):

  """ @brief  Attempts to compute the order r of g mod N, or a positive integer
              multiple thereof, given a frequency j yielded by the quantum part
              of Shor's order-finding algorithm, by using the post-processing
              algorithms described in detail in [E22p].

      @remark   This convenience function simply calls solve_j_for_r() with g
                setup by calling IntegerModRingMulSubgroupElement(g, N).

      [E22p] Ekerå, M.: "On the success probability of quantum order finding".
                        ArXiv 2201.07791 (2022).

      The idea is to try to solve not only j, but also j ± 1, .., j ± B, for r,
      with the aim of solving an optimal frequency j0(z) for r, for z the peak
      index on [0, r). Provided

        - that j0(z) is solved for r,

        - that d = gcd(r, z) is cm-smooth
          (for the definition of cm-smooth in [E22p]),

        - that l is selected as required by the solution method
          (see the documentation for parameter l below), and

        - that the accept_multiple flag is set to False,

      the order r will be found by this function.

      Note that this function does not implement meet-in-the-middle-techniques,
      although it is noted in [E23p] that it is possible to use such techniques
      to speed up the post-processing also for order finding.

      [E23p] Ekerå, M.: "On the success probability of the quantum algorithm for
                         the short DLP". ArXiv 2309.01754 (2023).

      @param j  The frequency j yielded by the quantum order-finding algorithm.

      @param m  A positive integer m such that r < 2^m.

      @param l  A positive integer l <= m, such that m+l is the length of the
                control register in the quantum order-finding algorithm.

                If method is set to SolutionMethods.CONTINUED_FRACTIONS_BASED or
                SolutionMethods.LATTICE_BASED_SHORTEST_VECTOR, it is required
                that r^2 < 2^(m+l) or else r may not be found.

                If method is set to SolutionMethods.LATTICE_BASED_ENUMERATE, it
                is possible to select l = m - Delta for some Delta in [0, m), at
                the expense of enumerating at most 6 * sqrt(3) * (2 ** Delta)
                lattice vectors for each offset in j considered.

      @param g  The group element g of order r modulo N.

      @param N  The modulus N.

      @param c  A parameter c >= 1 that specifies the maximum size of the
                missing cm-smooth component d in r = d * r_tilde when solving j
                for r, for cm-smooth as defined in [E22p].

                As is explained in [E22p], increasing c increases the success
                probability, at the expense of increasing the runtime.

      @param B  A bound B >= 0 on the offset in j. If B > 0, this function tries
                to solve not only j, but also j ± 1, .., j ± B, for r, or for a
                positive integer multiple of r.


      @param accept_multiple  A flag that may be set to True to indicate that
                              only a positive integer multiple of r is sought.
                              If set to True, this function returns as soon as
                              it finds r such that g^r = 1.

      @param method   An enumeration entry from the SolutionMethods class that
                      specifies the method to use to solve j for r. For further
                      details, see the documentation for the SolutionMethods
                      class.

      @param verbose  A flag that may be set to True to print intermediary
                      results and status updates when executing the
                      post-processing algorithm.

      @param opt_speculative  A flag that may be set to True to indicate that
                              Algorithm 2 in [E22p] should be used instead of
                              Algorithm 3 to find the missing cm-smooth
                              component of r. In most cases, Algorithm 2 is
                              faster than Algorithm 3, but in the worst case
                              Algorithm 2 is a lot slower than Algorithm 3. For
                              further details, see [E22p].

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

      @return   If the accept_multiple flag is set to False, the order r of g
                is returned with probability >= P, for P as given by the lower
                bound in [E22p] (provided that the opt_isolate_peak flag is set
                to False). Otherwise, None, or exceptionally a positive integer
                multiple of the order r, is returned. If the accept_multiple
                flag is set to True, some positive integer multiple of r is
                returned with probabilty >= P. Otherwise, None is returned. """

  g = IntegerModRingMulSubgroupElement(g, N);

  return solve_j_for_r(
    j = j,
    m = m,
    l = l,
    g = g,
    c = c,
    B = B,
    accept_multiple = accept_multiple,
    method = method,
    verbose = verbose,
    opt_speculative = opt_speculative,
    opt_isolate_peak = opt_isolate_peak);