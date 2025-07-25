""" @brief  A module for enumerating vectors in lattices. """

from .gram_schmidt import gram_schmidt;

from ..matrices import dimensions;
from ..matrices import solve_left;

from ..vectors import norm2;

from ...utils.timeout import Timeout;

from typing import Generator;

from gmpy2 import mpz;
from gmpy2 import mpq;
from gmpy2 import mpfr;

def enumerate(
  B : list[list[int | mpz]],
  R2 : int | mpz | mpq | mpfr,
  t : list[int | mpz | mpq | float | mpfr] | None = None,
  timeout : int | None | Timeout = None,
  gs : list[list[list[int | mpz | mpq | mpfr]]] | None = None,
  precision : int | None = None) -> \
    Generator[list[int | mpz], None, None]:

  """ @brief  Enumerates all vectors in the lattice L generated by the rows
              of the basis B that are within distance R of the vector t.

      The n x d integer basis matrix B is represented as as list of lists
      [b_1, ..., b_n], for b_i = [b_i1, ..., b_id] for i = 1, ..., n the n row
      vectors that make up B where b_ij is of type int or mpz.

      This function currently requires B to be square and to have full row rank
      since it calls solve_left(B, t) which takes the inverse of B to solve.
      These requirement may be relaxed in the future.

      @param B  The n x d basis matrix B = [b_1, ..., b_n], where
                b_i = [b_i1, ..., b_id] for i = 1, ..., n represent the n row
                vectors that make up the basis.

      @param R2  The square of the enumeration radius R.

      @param t  The vector t = [t_1, ..., t_d] in span(L). May be set to None,
                in which case t = [0, ..., 0].

      @param timeout  A timeout in seconds after which a TimeoutError will be
                      raised and the enumeration aborted. The timeout is
                      represented as integer or as an instance of Timeout. May
                      be set to None, as is the default, in which case no
                      timeout is enforced.

                      Note this function incrementally yields the vectors that
                      are found during the enumeration. The timeout specified
                      here is respect to the time elapsed since enumerate() was
                      first called. It hence includes any time spent processing
                      the results yielded by this function until control returns
                      to this function.

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
                        case the projection factors are represented as exact
                        quotients.

                        Note that this parameter only has an effect if gs is
                        set to None as the Gram–Schmidt orthogonalized basis Bs
                        and the associated matrix M of Gram–Schmidt projection
                        factors are otherwise pre-computed.

      return  A list [x_1, ..., x_m] of the m vectors found, where each vector
              x_i = [x_i1, ..., x_id] for 1, ..., m. This function yields one
              vector at a time allowing you to enumerate over the vectors. """

  # Initial setup.
  timeout = Timeout.parse(timeout);
  timeout.check();

  (n, d) = dimensions(B);

  if None == gs:
    [Bs, M] = gram_schmidt(B, precision = precision);
  else:
    [Bs, M] = gs;

  timeout.check();

  # 1
  if None == t:
    t = [0 for _ in range(d)];
    ct = [0 for _ in range(n)];
  else:
    ct = solve_left(B, t);

  cx = [0 for _ in range(n)];

  # 2
  def enumerate_inner(cx, ct, k):

    # Check the timeout.
    timeout.check();

    # 2.1
    if k == 0:

      # 2.1.1
      x = [sum([cx[i] * B[i][j] for i in range(n)]) for j in range(d)];
      yield x;

      return;

    # 2.2
    cx[k - 1] = 0;

    # 2.3
    a = ct[k - 1] - sum([M[i - 1][k - 1] * (cx[i - 1] - ct[i - 1]) \
                         for i in range(k + 1, n + 1)]);

    pi_k_1 = 0;
    for j in range(k + 1, n + 1):
      term = ((cx[j - 1] - ct[j - 1]) + \
                sum([M[i - 1][j - 1] * (cx[i - 1] - ct[i - 1]) \
                       for i in range(j + 1, n + 1)])) ** 2;
      term *= norm2(Bs[j - 1]);
      pi_k_1 += term;

    B2 = (R2 - pi_k_1) / norm2(Bs[k - 1]);
    if B2 < 0:
      return; # Return nothing.

    cx[k - 1] = round(a)
    tmp = abs(cx[k - 1] - a);
    if tmp * tmp <= B2:
      yield from enumerate_inner(cx.copy(), ct.copy(), k - 1);

    o = 0;
    proceed = True;

    while proceed:

      # Check the timeout.
      timeout.check();

      o = o + 1;
      proceed = False;

      cx[k - 1] = round(a) + o;
      tmp = abs(cx[k - 1] - a);
      if tmp * tmp <= B2:
        yield from enumerate_inner(cx.copy(), ct.copy(), k - 1);
        proceed = True;

      cx[k - 1] = round(a) - o;
      tmp = abs(cx[k - 1] - a);
      if tmp * tmp <= B2:
        yield from enumerate_inner(cx.copy(), ct.copy(), k - 1);
        proceed = True;

  # 3
  yield from enumerate_inner(cx, ct, n);
