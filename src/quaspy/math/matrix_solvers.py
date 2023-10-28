""" @brief  A module for solving c B = o for c, given a 2 x 2 integer matrix B,
            and an integer row vector o such that c B = o. """

def solve_left(B, o):

  """ @brief  Given a 2 x 2 integer matrix B, and an integer row vector o such
              that c B = o, this function returns the integer row vector c, or
              None if the equation c B = o has no integer solution

      @param B  The 2 x 2 integer matrix B = [b1, b2], where b1 = [b_11, b_22]
                and b2 = [b_21, b_22] are the row vectors that make up B.

      @param o  The integer row vector o = [o1, o2].

      @return   The integer row vector c = [c1, c2] such that c B = o, or None
                if the equation c B = o has no integer solution. """

  # Compute the determinant of B.
  det_B = B[0][0] * B[1][1] - B[0][1] * B[1][0];
  if 0 == det_B:
    # The matrix B does not have an inverse.
    return None;

  # Compute nu0.
  tmp = B[1][1] * o[0] - B[1][0] * o[1];
  if (tmp % det_B) != 0:
    # The equation c B = o has no integer solution
    return None;
  nu0 = tmp // det_B;

  # Compute nu1.
  tmp = -B[0][1] * o[0] + B[0][0] * o[1];
  if (tmp % det_B) != 0:
    # The equation c B = o has no integer solution
    return None;
  nu1 = tmp // det_B;

  # Return result.
  return [nu0, nu1];