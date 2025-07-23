""" @brief  A module for Gram–Schmidt orthogonalization. """

import gmpy2;

from gmpy2 import mpz;
from gmpy2 import mpq;
from gmpy2 import mpfr;

from ..matrices import dimensions;

from copy import deepcopy;

def gram_schmidt(
  B : list[list[int | mpz]],
  precision : int | None = None) -> \
    list[list[list[int | mpz | mpq | mpfr]]]:

  """ @brief  Returns the Gram–Schmidt orthogonalization Bs of an n x d basis
              matrix B, and the associated matrix M of Gram–Schmidt projection
              factors which is such that B = M Bs.

      The integer matrix B = [b_1, ..., b_n], where b_i = [b_i1, ..., b_id] for
      i = 1, ..., n are the n row vectors of B. The matrices Bs and M are
      represented in analogy with how B is represented.

      @param B  The matrix B.

      @param precision  The precision to use when computing the Gram–Schmidt
                        projection factors in M. May be set to None, as is the
                        default, in which case the projection factors are
                        represented as exact quotients.

      @return The pair [Bs, M], where Bs is the Gram–Schmidt orthogonalization
              of B and M is the matrix of Gram–Schmidt projection factors. """

  (n, d) = dimensions(B);

  B = [[mpz(x) for x in row] for row in B];

  Bs = deepcopy(B);

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

    for i in range(2, n+1):
      for j in range(i-1, 0, -1):
        M[i-1][j-1] = compute_mu(i-1, j-1);

        for k in range(d):
          Bs[i-1][k] = Bs[i-1][k] - M[i-1][j-1] * Bs[j-1][k];

    return [Bs, M];