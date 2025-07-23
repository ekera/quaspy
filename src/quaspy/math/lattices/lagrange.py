""" @brief  A module for computing Lagrange-reduced bases for two-dimensional
            lattices. """

from gmpy2 import mpz;
from gmpy2 import mpq;

from gmpy2 import t_divmod;

from ..vectors import norm2;

def lagrange(
    A : list[list[int | mpz]],
    multiples = None) -> list[list[list[int | mpz]]]:

  """ @brief  Returns the Lagrange-reduced basis for a 2 x 2 basis matrix A.

      The basis matrix A is represented as as list [u, v], where u = [u_1, u_2]
      and v = [v_1, v_2] represent the two row vectors that make up the basis.

      @param A  The 2 x 2 basis basis matrix A = [u, v], where u = [u_1, u_2]
                and v = [v_1, v_2] represent the two row vectors of the basis.

      @param multiples  Row multiples of the form [[m_uu, m_uv], [m_vu, m_vv]]
                        such that [m_uu * u + m_uv * v, m_vu * u + m_vv * v] is
                        a close to Lagrange-reduced basis for A, or None as is
                        the default if no such side information is available.

                        If row multiples are provided, they are used to start
                        the reduction process with a close to reduced basis,
                        typically speeding up the reduction process.

                        The matrix representing the row multiple must have full
                        rank if row multiples are provided. For as long as this
                        basic requirement is met, a Lagrange-reduced basis for A
                        will be returned even if the row multiples are way off.

      @return   The pair [A', multiples'], where A' = [u', v'] is a
                Lagrange-reduced basis for A = [u, v], and multiples' is of the
                form [[m'_uu, m'_uv], [m'_vu, m'_vv]] and A' = [u', v'] =
                [m'_uu * u + m'uv * v, m'_vu * u + m'_vv * v]. """

  # Input.
  [u, v] = A;

  u = [mpz(u[0]), mpz(u[1])];
  v = [mpz(v[0]), mpz(v[1])];

  if None == multiples:
    u_multiples = [mpz(1), mpz(0)];
    v_multiples = [mpz(0), mpz(1)];
  else:
    # Setup the combination of multiples.
    u_multiples = [mpz(multiples[0][0]), mpz(multiples[0][1])];
    v_multiples = [mpz(multiples[1][0]), mpz(multiples[1][1])];

    # Apply this combination of multiples to the basis for A.
    tmp = [u_multiples[0] * u[0] + u_multiples[1] * v[0],
           u_multiples[0] * u[1] + u_multiples[1] * v[1]];

    v = [v_multiples[0] * u[0] + v_multiples[1] * v[0],
         v_multiples[0] * u[1] + v_multiples[1] * v[1]];

    u = tmp;

  u_norm2 = norm2(u);
  v_norm2 = norm2(v);

  if u_norm2 < v_norm2:
    tmp = u;
    u = v;
    v = tmp;

    tmp = u_norm2;
    u_norm2 = v_norm2;
    v_norm2 = tmp;

    tmp = u_multiples;
    u_multiples = v_multiples;
    v_multiples = tmp;

  # Actual algorithm.
  while v_norm2 <= u_norm2:
    # Project u onto v and round to the closest integer.
    tmp = u[0] * v[0] + u[1] * v[1];
    [q, r] = t_divmod(tmp, v_norm2);

    # Account for rounding.
    if tmp >= 0:
      if  2 * r >= v_norm2:
        q += 1; # Round.
    else:
      if -2 * r >= v_norm2:
        q -= 1; # Round.

    # Subtract this projected component from u to make u shorter; swap u and v.
    w = [u[0] - q * v[0], u[1] - q * v[1]];
    w_norm2 = norm2(w);
    u = v; u_norm2 = v_norm2;
    v = w; v_norm2 = w_norm2;

    # Update the multiples accordingly, and reflect the above swap.
    u_multiples[0] -= q * v_multiples[0];
    u_multiples[1] -= q * v_multiples[1];

    tmp = u_multiples;
    u_multiples = v_multiples;
    v_multiples = tmp;

  # Output.
  return [[u, v], [u_multiples, v_multiples]];


def is_lagrange_reduced(
    A : list[list[int | mpz]]) -> bool:

  """ @brief  Returns True if the 2 x 2 basis matrix A is Lagrange-reduced,
              False otherwise.

      The basis matrix A is represented as as list [a1, a2], where a1, a2 are
      row vectors such that a1 = [a_11, a_12] and a2 = [a_21, a_22].

      @param A  The 2 x 2 basis basis matrix A.

      @return   True if A is Lagrange-reduced, False otherwise. """

  [u, v] = A;

  u = [mpz(u[0]), mpz(u[1])];
  v = [mpz(v[0]), mpz(v[1])];

  # It must hold that |u| <= |v|.
  if norm2(u) > norm2(v):
    return False;

  # It must hold that |μ| <= 1/2.
  abs_mu = mpq(abs(u[0] * v[0] + u[1] * v[1]), norm2(u));

  return abs_mu <= 1/2;