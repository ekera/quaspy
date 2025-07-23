""" @brief  A module for sampling random integers. """

from secrets import randbelow;

def sample_integer(B : int) -> int:

  """ @brief  Returns an integer selected uniformly at random from [0, B).

      @remark   This function calls randbelow() to select an integer uniformly
                at random from [0, B). In practice, randbelow() returns an
                integer that may be conjectured to be indistinguishable from an
                integer that is selected uniformly at random from [0, B).

      @param B  The upper bound B.

      @return   An integer selected uniformly at random from [0, B). """

  if B < 1:
    raise Exception("Error: Incorrect parameters: Pick B >= 1.");

  return randbelow(B);


def sample_l_bit_integer(l : int) -> int:

  """ @brief  Returns an l-bit integer selected uniformly at random from the set
              of all such integers.

      @remark   This function calls sample_integer() to select an integer
                uniformly at random from [0, B) for B = 2^(l-1). In practice,
                sample_integer() returns an integer that may be conjectured to
                be indistinguishable from an integer that is selected uniformly
                at random from [0, B).

      @param l  The bit length l of the prime to sample.

      @return   An l-bit integer selected uniformly at random from the set of
                all such integers. """

  if l < 1:
    raise Exception("Error: Incorrect parameters: Pick l >= 1.");

  return 2 ** (l - 1) + sample_integer(2 ** (l - 1));