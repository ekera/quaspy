""" @brief  A module for Gram–Schmidt orthogonalization. """

from .projections import proj;

def gram_schmidt(A):

  """ @brief  Returns the Gram–Schmidt orthogonalization of a 2 x 2 matrix A.

      The basis matrix A is represented as a list [a1, a2], where
      a = [a_11, a_22] and a2 = [a_21, a_22] are the two row vectors that form
      the basis. It is required that A has integer entries.

      @param A  The 2 x 2 integer basis matrix A = [a1, a2], where
                a1 = [a_11, a_22] and a2 = [a_21, a_22] are the two row vectors
                that form the basis.

      @return The Gram–Schmidt orthogonalization of A. """

  # Extract vectors.
  [a1, a2] = A;

  # First vector.
  b1 = a1;

  # Second vector.
  mu = proj(b1, a2);

  b2 = [a2[0] - mu * b1[0], a2[1] - mu * b1[1]];

  # Return result.
  return [b1, b2];