""" @brief  A module for computing projections. """

from gmpy2 import mpq;

from .norms import norm2;

def proj(u, v):

  """ @brief  Returns the projection of the two-dimensional integer vector v
              onto the two-dimensional integer vector u.

      @param u  The two-dimensional integer vector u = [u1, u2].

      @param v  The two-dimensional integer vector v = [v1, v2].

      @return   The projection of the two-dimensional integer vector v onto the
                two-dimensional integer vector u. """

  return mpq(u[0] * v[0] + u[1] * v[1], norm2(u));