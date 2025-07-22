""" @brief  A module for matrices. """

from typing import Any;

from gmpy2 import mpz;
from gmpy2 import mpq;

from copy import deepcopy;


def dimensions(B : list[list[Any]]):

  """ @brief  Returns the dimensions of the matrix B.

      @param B  The matrix B.

      @return   The dimensions (n, d) of the n x d-dimensional matrix B. """

  n = len(B);
  d = len(B[0]);

  return (n, d);


def invert(B : list[list[int | mpz | mpq]]) -> list[list[int | mpz | mpq]]:

  """ @brief  Inverts a full-rank rational n x n-matrix B.

      @param B  The matrix B = [b_1, ..., b_n], for b_i = [b_i1, ..., b_in] for
                i = 1, ..., n the n rows of B. The entries b_ij must be of type
                int or mpz for integer B, or of type mpq for rational B.

      @return   The inverse of the matrix B. """

  (n, d) = dimensions(B);

  if n != d:
    raise ValueError("Error: The matrix B must be a square matrix.");

  # Duplicate B.
  B = deepcopy(B);

  # Set B^-1 to I.
  B_inv = [[0 for _ in range(d)] for _ in range(n)];
  for i in range(n):
    B_inv[i][i] = 1;

  # The idea is now to apply row manipulations on B to transform it to I, and to
  # apply the exact same row manipulations to B_inv to transform it to B^-1.

  for column in range(n):

    # Find the first unused row in B that is non-zero in the column.
    row = column;
    while row < n:
      if B[row][column] != 0:
        break;
      row += 1;

    if row == n:
      raise ValueError("Error: The matrix B is not full rank.");

    # Force the entry in said row and column to one.
    c = mpq(1, B[row][column]);
    for j in range(d):
      B[row][j] *= c;
      B_inv[row][j] *= c;

    # Force the entries in all other rows in said column to zero.
    for i in range(n):
      if row == i:
        continue;

      c = B[i][column];
      for j in range(d):
        B[i][j] -= c * B[row][j];
        B_inv[i][j] -= c * B_inv[row][j];

    # Move the row to the right position if necessary.
    if row != column:
      tmp = B[row];
      B[row] = B[column];
      B[column] = tmp;

      tmp = B_inv[row];
      B_inv[row] = B_inv[column];
      B_inv[column] = tmp;

  # Return the inverse.
  return B_inv;


def solve_left(
  B : list[list[int | mpz | mpq]],
  t : list[int | mpz | mpq],
  B_inv : list[list[int | mpz | mpq]] | None = None) -> \
    list[int | mpz | mpq]:

  """ @brief  Given a full-rank rational n x n matrix B, and an n-dimensional
              rational row vector t such that c B = t where c is also a rational
              n-dimensional row vector, this function returns c.

      @param B  The matrix B = [b_1, ..., b_n], for b_i = [b_i1, ..., b_in] for
                i = 1, ..., n the n rows of B. The entries b_ij for i, j in
                1, ..., n must be of type int, mpz or mpq.

      @param t  The n-dimensional rational row vector t = [t_1, ..., t_n]. The
                entries t_i for i = 1, ..., n must be of type int, mpz or mpq.

      @param B_inv  The inverse of the matrix B as computed by calling the
                    function inverse(B), or None in which case this function
                    will call said function to compute the inverse.

      @return   The n-dimensional row vector c = [c_1, ..., c_n] such that
                c B = t. """

  (n, d) = dimensions(B);

  # Check that B is a square matrix.
  if n != d:
    raise Exception("Error: The matrix B must be a square matrix.");

  # We seek to solve c B = t for t. We do so simply by multiplying to the right
  # with the inverse B^-1 of B, yielding c = c B B^-1 = t B^-1.
  if None == B_inv:
    # Compute the inverse B^-1 of B.
    B_inv = invert(B);

  # Multiply t to the right by the inverse to form c = t B^-1.
  c = [sum([t[i] * B_inv[i][j] for i in range(n)]) for j in range(d)]

  # Return c.
  return c;