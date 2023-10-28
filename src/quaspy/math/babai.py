""" @brief  A module for Babai's nearest plane algorithm. """

from gmpy2 import mpz;

from .gram_schmidt import gram_schmidt;

from .projections import proj;

def babai(A, v, B = None):

  """ @brief  Uses Babai's nearest plane algorithm [Babai86] to find the vector
              in the lattice generated L by A that is closest to the vector v.

      [Babai86] L. Babai. On Lovász’ lattice reduction and the nearest lattice
                          point problem. Combinatorica 6(1), pp. 1–13 (1986).

      The matrix A is represented as a list [a1, a2], where a = [a_11, a_22] and
      a2 = [a_21, a_22] are the two row vectors that form the basis. It is
      required that A is Lagrange reduced, and that A has integer entries.

      @param A  The 2 x 2 Lagrange-reduced integer basis matrix A = [a_1, a_2]
                for the lattice L, where a_1 = [a_11, a_22] and
                a_2 = [a_21, a_22] are the two row vectors that form the basis.

      @param v  The integer row vector v = [v_1, v_2].

      @param B  The Gram–Schmidt orthogonalization B of A, if available, or None
                in which case the Gram–Schmidt orthogonalization B is computed
                by this function. (If you plan on calling this function several
                times for the same lattice L, then you may wish to pre-compute B
                and pass B along in each call to this function.)

      @return   The vector in the lattice L spanned by A that is closest to v.
  """

  if None == B:
    # Compute the Gram–Schmidt reduced basis for A.
    B = gram_schmidt(A);

  # Extract vectors.
  [a1, a2] = A;
  [b1, b2] = B;

  w = v;

  # Second vector.
  c = proj(b2, w);
  c_round = mpz(round(c));
  u = [c_round * a2[0], c_round * a2[1]];
  w = [w[0] - (c - c_round) * b2[0] - c_round * a2[0], \
       w[1] - (c - c_round) * b2[1] - c_round * a2[1]];

  # First vector.
  c = proj(b1, w);
  c_round = mpz(round(c));
  u = [u[0] + c_round * a1[0], u[1] + c_round * a1[1]];
  # w = [w[0] - (c - c_round) * b1[0] - c_round * a1[0], \
  #      w[1] - (c - c_round) * b1[1] - c_round * a1[1]];

  # Return result.
  return u;
