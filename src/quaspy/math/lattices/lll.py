""" @brief  A module for computing delta-LLL-reduced bases. """

import gmpy2;

from gmpy2 import mpz;
from gmpy2 import mpq;
from gmpy2 import mpfr;

from .gram_schmidt import gram_schmidt;

from ..matrices import dimensions;
from ..vectors import norm2;

from ...utils.timeout import Timeout;

LLL_DEFAULT_DELTA = 0.99;

LLL_DEFAULT_REDUCED_PRECISION = 128;

def lll(
  A : list[list[int | mpz]],
  delta : float = LLL_DEFAULT_DELTA,
  timeout : int | None | Timeout = None,
  gs : bool = False,
  precision : int | None = None) -> list[list[int | mpz]]:

  """ @brief  Returns the delta-LLL-reduced basis for an n x d basis matrix A,
              where LLL is short for Lenstra–Lenstra–Lovász [LLL82].

      [LLL82] A. K. Lenstra, H. W. Lenstra Jr. and L. Lovász. Factoring
              polynomials with rational coefficients. Math. Ann. 26(4),
              pp. 515–534 (1982).

      The basis matrix A is represented as a list of lists [a_1, ..., a_n],
      where a_i = [a_i1, ..., a_id] for i = 1, ..., n are the n rows of A. It is
      required that A has integer entries of type int or mpz.

      @param A  The basis matrix A.

      @param delta  The delta parameter in the Lovász condition. Must be on
                    the interval (1/4, 1]. A polynomial runtime in the
                    dimension n of the lattice is only guaranteed for delta < 1.

      @param timeout  A timeout after which a TimeoutError will be raised and
                      the computation aborted.

                      The timeout may be represented as an integer specifying
                      the timeout in seconds, or as an instance of the Timeout
                      class. May be set to None, as is the default, in which
                      case no timeout is enforced.

      @param gs  A flag that may be set to True to return not only the
                 delta-LLL-reduced basis B of A, but also the Gram–Schmidt
                 orthogonalization Bs of B, and the matrix M of Gram–Schmidt
                 projection factors such that B = M Bs. Note that Bs and M are
                 always computed as a part of the LLL reduction process.

      @param precision  The precision to use when computing the Gram–Schmidt
                        projection factors in M. May be set to None, as is the
                        default, in which case the projection factors are
                        represented as exact quotients.

      @return   The delta-LLL-reduced basis B = [b_1, ..., b_n], where
                b_i = [b_i1, ..., b_id] for i = 1, ..., n represent the n row
                vectors that make up the basis, if gs is set to False.

                Otherwise, if gs is set to True, [B, [Bs, M]] is returned where
                B is as above, Bs is the Gram–Schmidt orthogonalization of B and
                M is the matrix of Gram–Schmidt projection factors. It holds
                that B = M Bs. The matrices Bs and M are represented in analogy
                with how A and B are represented. """

  # Initial setup.
  timeout = Timeout.parse(timeout);
  timeout.check();

  (n, d) = dimensions(A);

  B = [[mpz(x) for x in row] for row in A]; # = (b_1, ..., b_n)^T

  Bs = [[mpz(0) for _ in range(d)] for _ in range(n)]; # = (b_1^*, ..., b_n^*)^T

  M = [[mpz(0) for _ in range(d)] for _ in range(n)];
  for i in range(n):
    M[i][i] = mpz(1);

  with gmpy2.context(gmpy2.get_context()) as context:

    # Set the precision.
    if precision != None:
      context.precision = precision;

    def compute_mu(i, j):
      a = sum([B[i][k]  * Bs[j][k] for k in range(d)]);
      b = sum([Bs[j][k] * Bs[j][k] for k in range(d)]);

      if None == precision:
        return mpq(a, b);

      with gmpy2.context(gmpy2.get_context()) as context:

        # Set the precision.
        context.precision = precision;

        return mpfr(a) / mpfr(b);

    # 1
    Bs[0] = B[0].copy(); # b_1^* <- b_1.

    i = 2;

    # 2
    while i <= n:

      # Check the timeout.
      timeout.check();

      # 2.1: Perform size reduction for b_i and b_i^*:

      # 2.1.1
      Bs[i - 1] = [0 for _ in range(d)];

      # 2.1.2
      for j in range(i - 1, 0, -1):

        # Check the timeout.
        timeout.check();

        # 2.1.2.1
        M[i - 1][j - 1] = compute_mu(i - 1, j - 1);

        # 2.1.2.2
        if abs(M[i - 1][j - 1]) > 1 / 2:

          tmp = round(M[i - 1][j - 1]);

          # 2.1.2.2.1
          for k in range(d):
            B[i - 1][k] = B[i - 1][k] - tmp * B[j - 1][k];

          # 2.1.2.2.2
          M[i - 1][j - 1] = M[i - 1][j - 1] - tmp;

        # 2.1.2.3
        for k in range(d):
          Bs[i - 1][k] = Bs[i - 1][k] - M[i - 1][j - 1] * Bs[j - 1][k];

      # 2.1.3
      for k in range(d):
        Bs[i - 1][k] = Bs[i - 1][k] + B[i - 1][k];


      # 2.2: Enforce the Lovász condition for b_{i−1}^* and b_i^*:

      # 2.2.1
      if (delta - M[i - 1][i - 2] * M[i - 1][i - 2]) * norm2(Bs[i - 2]) <= \
        norm2(Bs[i - 1]):

        # 2.2.1.1
        i = i + 1;

      # 2.2.2
      else:

        # 2.2.2.1
        tmp = B[i - 2];
        B[i - 2] = B[i - 1];
        B[i - 1] = tmp;

        # 2.2.2.2
        if i == 2:

          # 2.2.2.2.1
          Bs[0] = B[0].copy();

        # 2.2.2.3
        else:

          # 2.2.2.3.1
          i = i - 1;

    if gs:
      return [B, [Bs, M]];

    return B;

def is_lll_reduced(
  B : list[list[int | mpz]],
  delta : float = LLL_DEFAULT_DELTA,
  gs : list[list[list[int | mpz | mpq | mpfr]]] | None = None,
  precision : int | None = None) -> bool:

  """ @brief  Checks if an n x d basis matrix B is delta-LLL-reduced, where LLL
              is short for Lenstra–Lenstra–Lovász [LLL82].

      [LLL82] A. K. Lenstra, H. W. Lenstra Jr. and L. Lovász. Factoring
              polynomials with rational coefficients. Math. Ann. 26(4),
              pp. 515–534 (1982).

      The basis matrix B is represented as a list of lists [b_1, ..., b_n],
      where b_i = [b_i1, ..., b_id] for i = 1, ..., n are the n rows of A. It is
      required that B has integer entries of type int or mpz.

      @param B    The basis matrix B.

      @param delta  The delta parameter in the Lovász condition. Must be on
                    the interval (1/4, 1].

      @param gs   The list [Bs, M], where Bs is the Gram–Schmidt orthogonalized
                  basis for B, and M is the associated matrix of Gram–Schmidt
                  projection factors, as returned by calling

                     [Bs, M] = gram_schmidt(B, precision = precision),

                  or None, in which case this function will make the above call.

                  If you plan to perform several enumerations in the same
                  lattice, then time may be saved by not re-computing Bs and M
                  for each call.

                  Note that the basis B is typically LLL-reduced by calling
                  lll() before calling this function. The lll() function can
                  return not only B but also Bs and M (since Bs and M are
                  incrementally computed as a part of the LLL reduction process)
                  allowing you to directly pass them along to this function.

      @param precision  The precision to use when computing the Gram–Schmidt
                        projection factors in M. May be set to None, in which
                        case projection factors are represented as exact
                        quotients.

                        Note that this parameter only has an effect if gs is
                        set to None as the Gram–Schmidt orthogonalized basis Bs
                        and the associated matrix M of Gram–Schmidt projection
                        factors are otherwise pre-computed.

      @return   True if the basis B is delta-LLL-reduced, False otherwise. """

  # Initial setup.
  (n, _) = dimensions(B);

  if None == gs:
    [Bs, M] = gram_schmidt(B, precision = precision);
  else:
    [Bs, M] = gs;

  # Step 1: Check that B is size reduced.
  for i in range(2, n + 1):
    for j in range(1, i):
      if abs(M[i - 1][j - 1]) > 1 / 2:
        return False;

  # Step 2: Check the Lovász condition.
  for i in range(1, n - 1):
    if delta * norm2(Bs[i - 1]) > \
         norm2(Bs[i]) + M[i][i - 1] ** 2 * norm2(Bs[i - 1]):
      return False;

  # Both conditions are met so the basis B is delta-LLL-reduced.
  return True;