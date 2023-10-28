""" @brief  A module for solving an optimal frequency j = j0(z) for
            r_tilde = r / d where d = gcd(r, z). """

import gmpy2;

from gmpy2 import mpfr;
from gmpy2 import mpz;

from gmpy2 import gcd;

from math import sqrt;
from math import floor;

from .collection import CandidateCollection;

from ......math.lagrange import lagrange;
from ......math.continued_fractions import continued_fractions;

from ......math.primes import prime_power_product;

from ......math.norms import norm2;

def solve_j_for_r_tilde_continued_fractions(j, m, l):

  """ @brief  For j = j0(z) an optimal frequency, for m such that r < 2^m and l
              such that r^2 < 2^(m+l), and any z in [0, r), this function
              recovers r_tilde = r / d where d = gcd(r, z) by expanding the
              quotient j / 2^(m+l) in continued fractions and returning the
              last denominator < 2^((m+l)/2) as described in [E22p].

      [E22p] Ekerå, M.: "On the success probability of quantum order finding".
                        ArXiv 2201.07791 (2022).

      By Lemma 6 in [E22p], this function is guaranteed to return r_tilde
      provided that the requirements on the input parameters are met.

      For further details, see Lemma 6, and Sect. 4 and App. B, of [E22p].

      @param j  An optimal frequency j0(z), for m and l as passed to this
                function, and for any z in [0, r).

                If some other frequency j on [0, 2^(m+l)) is passed to this
                function, it may or may not return r_tilde.

      @param m  A positive integer m such that r < 2^m.

      @param l  A positive integer l <= m, such that r^2 < 2^(m+l).

      @return   The last denominator < 2^((m+l)/2) in the continued fraction
                expansion of j / 2^(m+l). This denominator is guaranteed to be
                equal to r_tilde, provided that the requirements on the input
                parameters are met. """

  bound = mpfr(2) ** mpfr((m + l) / 2);

  r_tilde_candidates = continued_fractions(j, m, l, denominator_bound = bound);

  return r_tilde_candidates[-1];


def solve_j_for_r_tilde_lattice_svp(j, m, l, multiples = None):

  """ @brief  For j = j0(z) an optimal frequency, for m such that r < 2^m and l
              such that r^2 < 2^(m+l), and any z in [0, r), this function
              recovers r_tilde = r / d where d = gcd(r, z) by finding the
              shortest non-zero vector in a two-dimensional lattice L as
              described in [E22p].

      [E22p] Ekerå, M.: "On the success probability of quantum order finding".
                        ArXiv 2201.07791 (2022).

      More specifically, this function uses Lagrange's lattice basis reduction
      algorithm to find the shortest non-zero vector

        u = (rj - 2^(m+l) z, r / 2) / d

      in the lattice L spanned by (j, 1/2) and (2^(m+l), 0), and hence r_tilde,
      as the second component is r_tilde / 2. This function return r_tilde.

      By Lemma 7 in [E22p], provided that r < 2^m and r^2 < 2^(m+l), the second
      component of the shortest non-zero vector in L has r_tilde / 2 in its
      second component, up to sign of course.

      For further details, see Lemma 7, and Sect. 4 and App. C, of [E22p].

      @param j  An optimal frequency j0(z), for m and l as passed to this
                function, and for any z in [0, r).

                If some other frequency j on [0, 2^(m+l)) is passed to this
                function, it will also return an integer, that may or may not
                be equal to r_tilde.

      @param m  A positive integer m such that r < 2^m.

      @param l  A positive integer l <= m, such that r^2 < 2^(m+l).

      @param multiples  Row multiples to use as a starting point when computing
                        the Lagrange-reduced basis for the slightly scaled basis
                        A = [[j, 1], [2^(m+l+1), 0]], see the lagrange()
                        function for further details, or None, as is the
                        default, if no such multiples are available.

                        The idea is that when trying to solve not only j but
                        also j ± 1, .., j ± B for r_tilde, in the hope that this
                        will lead to an optimal frequency j0(z) being solved for
                        r_tilde, the row multiples that yield a Lagrange-reduced
                        basis for adjecent offsets in j are likely to be close.

                        For this reason, this function accepts row multiples as
                        input, and it furthermore returns as output the row
                        multiples that yield a Lagrange-reduced basis for the
                        value of j input, so that these can be fed back to this
                        function when processing j + 1 or j - 1, recursively.

      @return   The tuple [r_tilde_candidate, multiple'], where multiple' are
                the row multiples that yield a reduced basis for the value of
                j input, and r_tilde_candidate is equal to r_tilde provided that
                the requirements on the input parameters are met. """

  # Setup the basis matrix for the lattice L (scaled by a factor of two).
  A = [[mpz(2 * j), mpz(1)], [2 * (2 ** (m + l)), mpz(0)]];
  # Without scaling: A = [[mpz(offset_j), mpq(1, 2)], [pow2ml, mpz(0)]];

  # Reduce the basis matrix.
  [A, multiples] = lagrange(A, multiples = multiples);

  # Extract r_tilde from the reduced basis (scaled by a factor of two).
  r_tilde_candidate = abs(mpz(A[0][1]));
  # Without scaling: r_tilde_candidate = abs(mpz(2 * A[0][1]));

  return [r_tilde_candidate, multiples];


def solve_j_for_r_tilde_lattice_enumerate(
  j,
  m,
  l,
  g,
  g_pow_e_context = None,
  c = 1,
  accept_multiple = False,
  partial_exponentiation = True,
  filtered_r_tilde_candidates = None,
  dismissed_reduced_r_tilde_candidates = None,
  mu = 0,
  multiples = None,
  verbose = False):

  """ @brief  For j = j0(z) an optimal frequency, for m such that r < 2^m and
              l = m - Delta for some Delta in [0, m), and z in [0, r), this
              function recovers r_tilde = r / d where d = gcd(r, z), provided
              that d is cm-smooth, by enumerating at most 6 * sqrt(3) * 2^Delta)
              vectors in a two-dimensional lattice L as described in [E22p].

      [E22p] Ekerå, M.: "On the success probability of quantum order finding".
                        ArXiv 2201.07791 (2022).

      More specifically, this function uses Lagrange's lattice basis reduction
      algorithm to find a reduced basis for the lattice L spanned by (j, 1/2)
      and (2^(m+l), 0). It then enumerates all vectors within a ball of radius
      2^(m-1/2) in L centered at the origin to find

        u = (rj - 2^(m+l) z, r / 2) / d,

      and hence r_tilde, as the second component is r_tilde / 2.

      By Lemma 8 in [E22p], provided that r < 2^m and l = m - Delta, at most
      6 * sqrt(3) * 2^Delta vectors must be enumerated in L to find u and hence
      r_tilde, so if Delta is small then this method is efficient.

      In practice, as mentioned in [E22p], the leading constant in the above
      bound is not tight, and the enumeration can be optimized. Some of these
      optimizations have been integrated into this implementation, so the
      enumeration typically considers fewer points than the bound indicates.

      @remark   Unlike the solve_j_for_r_tilde_lattice_svp() function, and the
                solve_j_for_r_tilde_continued_fractions(), this function does
                check that the candidates x for r_tilde returned fulfill the
                requirement that e(x) * x is a positive multiple of r, where
                e(x) is cm-smooth by the definition of cm-smooth in [E22p].

                This is necessary as the enumeration may generate a
                comparatively, and it can also be done efficiently.

                A side effect of this test being performed is that this function
                will typically not find r_tilde = r / d if d is not cm-smooth.

      @remark   By default, this implementation uses an instance of the
                CandidateCollection to keep track of the collection of
                candidates for r_tilde found during the enumeration.

                As is explained in the documentation for the CandidateCollection
                class, this class keeps the collection reduced. In turn, this
                may lead to r_tilde not being returned as one of the candidates
                for r_tilde that represent the collection. This being said, one
                of the candidates x for r_tilde that represent the collection
                will then be such that r = e(x) * x, where e(x) is cm-smooth,
                and from such x it is always possible to correctly recover r.

                As is also explained below, you may input set() instead of
                CandidateCollection() to the filtered_r_tilde_candidates
                parameter to ensure all candidates for r_tilde are returned.

      For further details, see Lemma 8, and Sect. 4 and App. C, of [E22p].

      @param j  An optimal frequency j0(z), for m and l as passed to this
                function, and for z in [0, r).

      @param m  A positive integer m such that r < 2^m.

      @param l  A positive integer l <= m, such that m+l is the length of the
                control register in the quantum order-finding algorithm.

                As is explained in [E22p], it is possible to select l = m-Delta,
                for some Delta in [0, m), at the expense of enumerating at most
                6 * sqrt(3) * 2^Delta lattice vectors.

      @param g  The group element g of order r.

      @param g_pow_e_context  The pair [x, e] for x = g^e and e a maximal
                              cm-smooth product, or None, in which case this
                              function will compute x and e.

      @param c  A parameter c >= 1 that specifies the maximum size of the
                missing cm-smooth component d in r = d * r_tilde. As is
                explained in [E22p], increasing c increases the success
                probability, at the expense of increasing the runtime.

      @param accept_multiple  A flag that may be set to True to indicate that
                              only a positive integer multiple of r is sought.
                              If set to True, this function returns as soon as
                              it finds a candidate x for r_tilde such that
                              g^(d * x) = 1, for d a cm-smooth prime power
                              product.

      @param partial_exponentiation   A flag that may be set to True to indicate
                                      that x^s1[1] and x^s1[2] shall first be
                                      computed, and that the candidate for
                                      r_tilde for element i1 * s[1] + i2 * s[2]
                                      in the lattice shall then be computed by
                                      as (x1^i1) * (x2^i2).

      @param filtered_r_tilde_candidates  A collection of candidates x for
                                          r_tilde, such that g^(e(x) * x) = 1,
                                          for e(x) some cm-smooth prime power
                                          product.

                                          When e.g. solving not only j, but also
                                          j ± 1, .., j ± B, for r, for B some
                                          bound, then this parameter allows for
                                          accounting for candidates of r_tilde
                                          found in enumerations for previous
                                          offsets when enumerating for a given
                                          offset. In particular, it removes the
                                          need to perform exponentiations to
                                          test candidates x for r_tilde that
                                          have been found before, and checked
                                          to meet the requirement that
                                          g^(e(x) * x) = 1, for e(x) some
                                          cm-smooth prime power product.

                                          May be set to CandidateCollection(),
                                          as is the default when passing None,
                                          to start with an empty collection of
                                          candidates for r_tilde that is kept
                                          reduced, or to set(), to keep all
                                          candidates for r_tilde found during
                                          the enumeration in the set.

                                          The set of candidates for r_tilde
                                          returned by this function will be of
                                          the same type as the set input via
                                          this parameter.

      @param dismissed_reduced_r_tilde_candidates   A set of reduced candidates
                                                    for r_tilde that have
                                                    already been dismissed and
                                                    may be dismissed immediately
                                                    if encountered again.

                                                    May be set to set(), as is
                                                    the default when passing
                                                    None, if no candidates for
                                                    r_tilde have been dismissed.

      @param mu   A non-negative multiple of the order r used to reduce the
                  candidates for r_tilde found during the enumeration. May be
                  set to zero, as is the default, if no non-trivial multiple of
                  r is known.

      @param multiples  Row multiples to use as a starting point when computing
                        the Lagrange-reduced basis for the slightly scaled basis
                        A = [[j, 1], [2^(m+l+1), 0]], see the lagrange()
                        function for further details, or None, as is the
                        default, if no such multiples are available.

                        The idea is that when trying to solve not only j but
                        also j ± 1, .., j ± B for r_tilde, in the hope that this
                        will lead to an optimal frequency j0(z) being solved for
                        r_tilde, the row multiples that yield a Lagrange-reduced
                        basis for adjecent offsets in j are likely to be close.

                        For this reason, this function accepts row multiples as
                        input, and it furthermore returns as output the row
                        multiples that yield a Lagrange-reduced basis for the
                        value of j input, so that these can be fed back to this
                        function when processing j + 1 or j - 1, recursively.

      @param verbose  A flag that may be set to True to print intermediary
                      results and status updates when performing th
                      enumeration.

      @return   The tuple [filtered_r_tilde_candidates', [success,
                dismissed_reduced_r_tilde_candidates', mu', multiples']], where

                - filtered_r_tilde_candidates' is filtered_r_tilde_candidates
                  as input to this function updated with the additional
                  candidates for r_tilde found when enumerating the lattice
                  (these candidates have all passed the filtration step, and are
                   such that d * t is a multiple of r, where d is cm-smooth),

                - success is set to True if at least one candidate for r_tilde
                  passed the filter when enumerating the lattice, even if this
                  candidate was already in filtered_r_tilde_candidates, and
                  to False otherwise,

                - dismissed_reduced_r_tilde_candidates' is
                  dismissed_reduced_r_tilde_candidates as input to this function
                  updated with the additional reduced candidates for r_tilde
                  dismissed when enumerating the lattice,

                - mu' is an updated non-negative multiple of r, obtained by
                  letting mu = mu' and then updating mu' recursively by letting
                  mu' = gcd(mu', filted_r_tilde_candidate * e), where e is a
                  maximal cm-smooth product and x runs over all candidates for
                  r_tilde that pass the filtration step when enumerating the
                  lattice, and

                - multiple' are the row multiples that yield a reduced basis for
                  the value of j input. """


  # Enforce default. Note: We cannot use CandidateCollection() or set() as
  # defaults in the prototype because Python3 caches these collections/sets.
  if filtered_r_tilde_candidates == None:
    filtered_r_tilde_candidates = CandidateCollection();

  if dismissed_reduced_r_tilde_candidates == None:
    dismissed_reduced_r_tilde_candidates = set();

  # Sanity checks.
  if (m <= 0) or (l < 0) or (l > m):
    raise Exception("Error: Incorrect parameters");

  # Setup precision.
  swapped_out_precision = gmpy2.get_context().precision;
  gmpy2.get_context().precision = 53;

  # Compute Delta given m and l.
  Delta = m - l;

  # Setup the basis matrix for the lattice L (scaled by a factor of two).
  A = [[mpz(2 * j), mpz(1)], [2 * mpz(2 ** (m + l)), mpz(0)]];
  # Without scaling: A = [[mpz(j), mpq(1, 2)], [mpz(2 ** (m + l), mpz(0)]];

  # Reduce the basis matrix.
  [A, multiples] = lagrange(A, multiples = multiples);

  # Extract the shortest non-zero vector, denoted s1, and the shortest
  # non-zero vector that is linerly independent to s1, denoted s2.
  [s1, s2] = A;

  # Compute float representation of these vectors, since they may be large.
  s1f = [mpfr(x) for x in s1];
  s2f = [mpfr(x) for x in s2];

  # Compute the Gram-Schmidt-coefficient mu21, such that
  #   mu12 * s1 = component of s2 that is parallel to s1, and
  #   s2 - mu12 * s1 = component of s2 that is orthogonal to s1.
  mu12 = (s1f[0] * s2f[0] + s1f[1] * s2f[1]) / norm2(s1f);

  # Compute the parallel and orthogonal components of s2.
  s2f_parallel = [mu12 * s1f[0], mu12 * s1f[1]];
  s2f_orthogonal = [s2f[0] - s2f_parallel[0], s2f[1] - s2f_parallel[1]];

  if None == g_pow_e_context:
    # Form e.
    e = prime_power_product(c * m); # TBD: Pass e to this function alongside x.

    # Exponentiate g to e to form x.
    x = g ** e;
  else:
    [x, e] = g_pow_e_context;

    e = mpz(e);

  # The radius of the circle to enumerate. In [E22p], the radius of the circle
  # to enumerate is 2^(m - 1/2), which would imply radius2 = 2^(2m - 1). This
  # bound stems from the fact that the target vector is
  #
  #   | [alpha_0(z) / d, (r / 2) / d] | <= | [r/2, r/2] | =
  #     sqrt((r/2)^2 + (r/2)^2) = sqrt(r^2 / 2) = r / sqrt(2) <
  #       2^m / sqrt(2) = 2^(m-1/2),
  #
  # since d = gcd(z, r) >= 1, |alpha_0(z)| <= r/2 and r < 2^m.
  #
  # We have scaled the lattice by a factor of two, yielding instead
  #
  #   | [2 alpha_0(z) / d, r / d] | <= | [r, r] | =
  #     sqrt(r^2 + r^2) = sqrt(2 r^2) = sqrt(2) r <
  #       sqrt(2) 2^m = 2^(m+1/2),
  #
  # and (2^(m+1/2))^2 = 2^(2m+1) as below.
  radius2 = mpfr(mpz(2 ** (2 * m + 1)));

  # A flag that is set to True if a candidate was found and to False otherwise.
  success = False;

  # Setup mu.
  mu = mpz(mu);

  if norm2(s2f_orthogonal) >= radius2:
    # As is stated in [E22p], if | s2_orthogonal | >= radius^2, we have that
    # the second component of the shortest non-zero vector must be r_tilde / 2
    # (and we have scaled by factor of two, so the component is now r_tilde).

    r_tilde_candidate = abs(s1[1]);

    if r_tilde_candidate in filtered_r_tilde_candidates:
      success = True;
    else:
      # Use that mu is an r-multiple to reduce the candidate for r_tilde.
      reduced_r_tilde_candidate = gcd(r_tilde_candidate, mu);

      if (reduced_r_tilde_candidate in dismissed_reduced_r_tilde_candidates):
        # Dismiss the reduced candidate.
        if verbose:
          print("Dismissing:", r_tilde_candidate);
      else:
        # The reduced candidate has not already been dismissed.
        if verbose:
          print("Testing the reduced candidate:", \
            reduced_r_tilde_candidate);

        # Test the reduced candidate.
        if (x ** reduced_r_tilde_candidate) == 1:
          success = True;

          # Add r_tilde_candidate to the filtered candidates for r_tilde.
          filtered_r_tilde_candidates.add(r_tilde_candidate);

          # We know that reduced_r_tilde_candidate * e is a multiple of r,
          # so we may update mu to account for this fact:
          mu = gcd(reduced_r_tilde_candidate * e, mu);
        else:
          # Add reduced_r_tilde_candidate to the dismissed reduced
          # candidates for r_tilde to avoid future exponentiations.
          dismissed_reduced_r_tilde_candidates.\
            add(reduced_r_tilde_candidate);

    gmpy2.get_context().precision = swapped_out_precision;
    return [filtered_r_tilde_candidates,
            [success,
             dismissed_reduced_r_tilde_candidates,
             mu,
             multiples]];

  # Compute an upper bound B on the number of points to enumerate.
  B = floor(6 * sqrt(3) * (2 ** Delta));

  # Pre-compute 2^m for later comparisons.
  pow2mf = mpfr(mpz(2 ** m));

  # Storage for x_basis = [x ** s1[1], x ** s2[1]] that is precomputed upon
  # demand if the partial_exponentiation flag is set to True.
  x_basis = None;

  count = 0;

  i2 = 0;

  while True:
    # Check the condition on the radius.
    u2_orthogonalf = [i2 * s2f_orthogonal[0], i2 * s2f_orthogonal[1]];
    if norm2(u2_orthogonalf) > radius2:
      break;

    # Form u2f.
    u2f = [i2 * s2f[0], i2 * s2f[1]];

    # Form i1hat and search i1 = i1hat, i1hat ± 1, i1hat ± 2, ..
    i1 = i1hat = round(-mu12 * i2);

    # Optimization: Jump ahead in i1.
    uf = [u2f[0] + i1 * s1f[0], u2f[1] + i1 * s1f[1]];

    if s1f[1] >= 0:
      if uf[1] <= -pow2mf:
        i1 += int(floor((-uf[1] - pow2mf) // s1f[1]));
    else:
      if uf[1] >= pow2mf:
        i1 += int(floor((uf[1] - pow2mf) // -s1f[1]));
    # End of optimization.

    while True:
      # Check the condition on the radius.
      uf = [u2f[0] + i1 * s1f[0], u2f[1] + i1 * s1f[1]];
      if norm2(uf) > radius2:
        break;

      # Update the count with an additional candidate point. We only search over
      # positive i2, but if [i1, i2] is a candidate, then of course so is
      # [-i1, -i2], so unless i2 = 0 we count the candidate twice.
      if i2 != 0:
        count += 2;
      else:
        count += 1;

      # Check the candidate.
      if (not (i1 == i2 == 0)) and (0 < abs(uf[0]) < pow2mf) \
                               and (0 < abs(uf[1]) < pow2mf):

        # Compute r_tilde_candidate.
        r_tilde_candidate = abs(i1 * s1[1] + i2 * s2[1]);

        if r_tilde_candidate in filtered_r_tilde_candidates:
          success = True;

          if accept_multiple:
            gmpy2.get_context().precision = swapped_out_precision;
            return [filtered_r_tilde_candidates,
                    [success,
                     dismissed_reduced_r_tilde_candidates,
                     mu,
                     multiples]];
        else:
          # Use that mu is an r-multiple to reduce the candidate for r_tilde.
          reduced_r_tilde_candidate = gcd(r_tilde_candidate, mu);

          if (reduced_r_tilde_candidate in \
            dismissed_reduced_r_tilde_candidates):
            # Dismiss the reduced candidate.
            if verbose:
              print("Dismissing:", r_tilde_candidate);
          else:
            # The reduced candidate has not already been dismissed.
            if verbose:
              print("Testing the candidate:", i1, i2, \
                reduced_r_tilde_candidate, r_tilde_candidate);

            # Test the reduced candidate.
            if partial_exponentiation:
              if x_basis == None:
                  x_basis = [x ** s1[1], x ** s2[1]];

              z = (x_basis[0] ** i1) * (x_basis[1] ** i2);
            else:
              z = x ** reduced_r_tilde_candidate;

            if z == 1:
              success = True;

              # Add r_tilde_candidate to the filtered candidates for r_tilde.
              filtered_r_tilde_candidates.add(r_tilde_candidate);

              if accept_multiple:
                gmpy2.get_context().precision = swapped_out_precision;
                return [filtered_r_tilde_candidates,
                        [success,
                         dismissed_reduced_r_tilde_candidates,
                         mu,
                         multiples]];

              # We know that reduced_r_tilde_candidate * e is a multiple of r,
              # so we may update mu to account for this fact:
              mu = gcd(reduced_r_tilde_candidate * e, mu);
            else:
              # Add reduced_r_tilde_candidate to the dismissed reduced
              # candidates for r_tilde to avoid future exponentiations.
              dismissed_reduced_r_tilde_candidates.\
                add(reduced_r_tilde_candidate);

      if s1f[0] >= 0:
        if uf[0] >=  pow2mf:
          break; # There is no point in continuing.
      else:
        if uf[0] <= -pow2mf:
          break; # There is no point in continuing.

      if s1f[1] >= 0:
        if uf[1] >=  pow2mf:
          break; # There is no point in continuing.
      else:
        if uf[1] <= -pow2mf:
          break; # There is no point in continuing.

      i1 += 1;

    i1 = i1hat - 1;

    # Optimization: Jump ahead in i1.
    uf = [u2f[0] + i1 * s1f[0], u2f[1] + i1 * s1f[1]];

    if s1f[1] >= 0:
      if uf[1] <= -pow2mf:
        i1 -= int(floor((-uf[1] - pow2mf) // s1f[1]));
    else:
      if uf[1] >= pow2mf:
        i1 -= int(floor((uf[1] - pow2mf) // -s1f[1]));
    # End of optimization.

    while True:
      # Check the condition on the radius.
      uf = [u2f[0] + i1 * s1f[0], u2f[1] + i1 * s1f[1]];
      if norm2(uf) > radius2:
        break;

      # Update the count with an additional candidate point. We only search
      # over positive i2, but if [i1, i2] is a candidate, then of course so
      # is [-i1, -i2], so unless i2 = 0 we count the candidate twice.
      if i2 != 0:
        count += 2;
      else:
        count += 1;

      # Check the candidate.
      if (not (i1 == i2 == 0)) and (0 < abs(uf[0]) < pow2mf) \
                               and (0 < abs(uf[1]) < pow2mf):

        # Compute r_tilde_candidate.
        r_tilde_candidate = abs(i1 * s1[1] + i2 * s2[1]);

        if r_tilde_candidate in filtered_r_tilde_candidates:
          success = True;

          if accept_multiple:
            gmpy2.get_context().precision = swapped_out_precision;
            return [filtered_r_tilde_candidates,
                    [success,
                     dismissed_reduced_r_tilde_candidates,
                     mu,
                     multiples]];
        else:
          # Use that mu is an r-multiple to reduce the candidate for r_tilde.
          reduced_r_tilde_candidate = gcd(r_tilde_candidate, mu);

          if (reduced_r_tilde_candidate in \
            dismissed_reduced_r_tilde_candidates):
            # Dismiss the reduced candidate.
            if verbose:
              print("Dismissing:", r_tilde_candidate);
          else:
            # The reduced candidate has not already been dismissed.
            if verbose:
              print("Testing the candidate:", i1, i2, \
                reduced_r_tilde_candidate, r_tilde_candidate);

            # Test the reduced candidate.
            if partial_exponentiation:
              if x_basis == None:
                  x_basis = [x ** s1[1], x ** s2[1]];

              z = (x_basis[0] ** i1) * (x_basis[1] ** i2);
            else:
              z = x ** reduced_r_tilde_candidate;

            if z == 1:
              success = True;

              # Add r_tilde_candidate to the filtered candidates for r_tilde.
              filtered_r_tilde_candidates.add(r_tilde_candidate);

              if accept_multiple:
                gmpy2.get_context().precision = swapped_out_precision;
                return [filtered_r_tilde_candidates,
                        [success,
                         dismissed_reduced_r_tilde_candidates,
                         mu,
                         multiples]];

              # We know that reduced_r_tilde_candidate * e is a multiple of r,
              # so we may update mu to account for this fact:
              mu = gcd(reduced_r_tilde_candidate * e, mu);
            else:
              # Add reduced_r_tilde_candidate to the dismissed reduced
              # candidates for r_tilde to avoid future exponentiations.
              dismissed_reduced_r_tilde_candidates.\
                add(reduced_r_tilde_candidate);

      if s1f[0] <= 0:
        if uf[0] >=  pow2mf:
          break; # There is no point in continuing.
      else:
        if uf[0] <= -pow2mf:
          break; # There is no point in continuing.

      if s1f[1] <= 0:
        if uf[1] >=  pow2mf:
          break; # There is no point in continuing.
      else:
        if uf[1] <= -pow2mf:
          break; # There is no point in continuing.

      i1 -= 1;

    # Try next i2.
    i2 += 1;

    # Sanity check.
    if count >= B:
      raise Exception("Error: Enumerated more vectors than expected.");

  gmpy2.get_context().precision = swapped_out_precision;
  return [filtered_r_tilde_candidates,
          [success,
           dismissed_reduced_r_tilde_candidates,
           mu,
           multiples]];